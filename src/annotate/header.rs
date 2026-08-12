use std::{
    collections::{HashMap, HashSet},
    fs,
    path::{Path, PathBuf},
};

use noodles_vcf::{
    Header,
    variant::{
        RecordBuf,
        record_buf::{Samples, samples::Keys},
    },
};
use rsomics_common::{Context, Result, RsomicsError};

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct HeaderOptions {
    pub(crate) appended: Vec<String>,
    pub(crate) remove: Option<String>,
    pub(crate) rename_chromosomes: Option<PathBuf>,
    pub(crate) rename_annotations: Option<PathBuf>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Removal {
    Id,
    Qual,
    Filters(Selection<String>),
    Info(Selection<String>),
    Format(Selection<String>),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Selection<T> {
    pub(crate) keep: bool,
    pub(crate) values: Vec<T>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct Renames {
    pub(crate) chromosomes: HashMap<String, String>,
    pub(crate) infos: HashMap<String, String>,
    pub(crate) formats: HashMap<String, String>,
    pub(crate) filters: HashMap<String, String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct HeaderPlan {
    appended: Vec<String>,
    removals: Vec<Removal>,
    renames: Renames,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
enum Namespace {
    Info,
    Format,
    Filter,
}

#[derive(Default)]
struct SelectionBuilder {
    mode: Option<bool>,
    all: bool,
    values: Vec<String>,
    seen: HashSet<String>,
}

impl HeaderPlan {
    pub(crate) fn bind(header: &Header, options: HeaderOptions) -> Result<Self> {
        let removals = parse_removals(options.remove.as_deref())?;
        let renames = Renames {
            chromosomes: read_chromosome_renames(options.rename_chromosomes.as_deref())?,
            ..read_annotation_renames(options.rename_annotations.as_deref())?
        };
        let plan = Self {
            appended: options.appended,
            removals,
            renames,
        };
        let mut candidate = header.clone();
        plan.prepare(&mut candidate)?;
        Ok(plan)
    }

    pub(crate) fn prepare(&self, header: &mut Header) -> Result<()> {
        self.prepare_with_retained_removals(header, false)
    }

    pub(crate) fn prepare_with_retained_removals(
        &self,
        header: &mut Header,
        retain: bool,
    ) -> Result<()> {
        let mut candidate = header.clone();
        append_header_lines(&mut candidate, &self.appended)?;
        validate_removals(&candidate, &self.removals)?;
        if !retain {
            remove_header_definitions(&mut candidate, &self.removals);
        }
        validate_renames(&candidate, &self.renames)?;
        rename_header_maps(&mut candidate, &self.renames);
        *header = candidate;
        Ok(())
    }

    pub(crate) fn apply(&self, record: &mut RecordBuf) -> Result<bool> {
        let mut changed = apply_removals(record, &self.removals);
        changed |= apply_renames(record, &self.renames)?;
        Ok(changed)
    }

    pub(crate) fn apply_renames(&self, record: &mut RecordBuf) -> Result<bool> {
        apply_renames(record, &self.renames)
    }
}

impl Selection<String> {
    fn removes(&self, value: &str) -> bool {
        let listed = self.values.iter().any(|candidate| candidate == value);
        if self.keep {
            !listed
        } else {
            self.values.is_empty() || listed
        }
    }
}

fn parse_removals(source: Option<&str>) -> Result<Vec<Removal>> {
    let Some(source) = source else {
        return Ok(Vec::new());
    };
    if source.trim().is_empty() {
        return Err(invalid("annotation removal list must not be empty"));
    }

    let mut id = false;
    let mut qual = false;
    let mut info = SelectionBuilder::default();
    let mut format = SelectionBuilder::default();
    let mut filter = SelectionBuilder::default();

    for raw in source.split(',').map(str::trim) {
        if raw.is_empty() {
            return Err(invalid("annotation removal list contains an empty item"));
        }
        let retained = raw.starts_with('^');
        let value = raw.strip_prefix('^').unwrap_or(raw);
        if value.eq_ignore_ascii_case("ID") {
            if retained || std::mem::replace(&mut id, true) {
                return Err(invalid(format!(
                    "invalid or duplicate annotation removal: {raw}"
                )));
            }
            continue;
        }
        if value.eq_ignore_ascii_case("QUAL") {
            if retained || std::mem::replace(&mut qual, true) {
                return Err(invalid(format!(
                    "invalid or duplicate annotation removal: {raw}"
                )));
            }
            continue;
        }
        if let Some(namespace) = bare_namespace(value) {
            if retained {
                return Err(invalid(format!(
                    "invalid annotation retained namespace: {raw}"
                )));
            }
            selection_builder(namespace, &mut info, &mut format, &mut filter).set_all(raw)?;
            continue;
        }
        let (namespace, tag) = typed_name(value)?;
        selection_builder(namespace, &mut info, &mut format, &mut filter)
            .insert(tag, retained, raw)?;
    }

    let mut removals = Vec::new();
    if id {
        removals.push(Removal::Id);
    }
    if qual {
        removals.push(Removal::Qual);
    }
    if let Some(selection) = info.finish() {
        removals.push(Removal::Info(selection));
    }
    if let Some(selection) = format.finish() {
        removals.push(Removal::Format(selection));
    }
    if let Some(selection) = filter.finish() {
        removals.push(Removal::Filters(selection));
    }
    Ok(removals)
}

impl SelectionBuilder {
    fn set_all(&mut self, raw: &str) -> Result<()> {
        if self.mode.is_some() || self.all {
            return Err(invalid(format!("conflicting annotation removal: {raw}")));
        }
        self.mode = Some(false);
        self.all = true;
        Ok(())
    }

    fn insert(&mut self, value: String, retained: bool, raw: &str) -> Result<()> {
        let retained = retained || self.mode == Some(true);
        if self.all || self.mode.is_some_and(|mode| mode != retained) {
            return Err(invalid(format!("conflicting annotation removal: {raw}")));
        }
        self.mode = Some(retained);
        if !self.seen.insert(value.clone()) {
            return Err(invalid(format!("duplicate annotation removal: {raw}")));
        }
        self.values.push(value);
        Ok(())
    }

    fn finish(self) -> Option<Selection<String>> {
        self.mode.map(|keep| Selection {
            keep,
            values: self.values,
        })
    }
}

fn selection_builder<'a>(
    namespace: Namespace,
    info: &'a mut SelectionBuilder,
    format: &'a mut SelectionBuilder,
    filter: &'a mut SelectionBuilder,
) -> &'a mut SelectionBuilder {
    match namespace {
        Namespace::Info => info,
        Namespace::Format => format,
        Namespace::Filter => filter,
    }
}

fn bare_namespace(value: &str) -> Option<Namespace> {
    if value.eq_ignore_ascii_case("INFO") || value.eq_ignore_ascii_case("INF") {
        Some(Namespace::Info)
    } else if value.eq_ignore_ascii_case("FORMAT") || value.eq_ignore_ascii_case("FMT") {
        Some(Namespace::Format)
    } else if value.eq_ignore_ascii_case("FILTER") {
        Some(Namespace::Filter)
    } else {
        None
    }
}

fn typed_name(value: &str) -> Result<(Namespace, String)> {
    let Some((prefix, tag)) = value.split_once('/') else {
        return Err(invalid(format!(
            "annotation removal requires INFO/, FORMAT/, or FILTER/: {value}"
        )));
    };
    let namespace = bare_namespace(prefix).ok_or_else(|| {
        invalid(format!(
            "unsupported annotation namespace in removal: {value}"
        ))
    })?;
    if !valid_tag(tag) {
        return Err(invalid(format!("invalid annotation tag: {value}")));
    }
    Ok((namespace, tag.to_owned()))
}

fn append_header_lines(header: &mut Header, lines: &[String]) -> Result<()> {
    for line in lines {
        if line.contains(['\n', '\r']) {
            return Err(invalid(
                "an appended VCF header definition spans multiple lines",
            ));
        }
        let file_format = header.file_format();
        let raw = format!(
            "##fileformat=VCFv{}.{}\n{line}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            file_format.major(),
            file_format.minor()
        );
        let parsed: Header = raw.parse().map_err(|error| {
            invalid(format!(
                "invalid appended VCF header definition {line:?}: {error}"
            ))
        })?;

        let populated = usize::from(!parsed.infos().is_empty())
            + usize::from(!parsed.formats().is_empty())
            + usize::from(!parsed.filters().is_empty())
            + usize::from(!parsed.contigs().is_empty());
        if populated != 1
            || parsed.infos().len() > 1
            || parsed.formats().len() > 1
            || parsed.filters().len() > 1
            || parsed.contigs().len() > 1
        {
            return Err(invalid(format!(
                "unsupported appended VCF header definition: {line}"
            )));
        }
        if let Some((name, value)) = parsed.infos().first() {
            if header.infos().contains_key(name) {
                return Err(invalid(format!(
                    "INFO/{name} already exists in the VCF header"
                )));
            }
            header.infos_mut().insert(name.clone(), value.clone());
        } else if let Some((name, value)) = parsed.formats().first() {
            if header.formats().contains_key(name) {
                return Err(invalid(format!(
                    "FORMAT/{name} already exists in the VCF header"
                )));
            }
            header.formats_mut().insert(name.clone(), value.clone());
        } else if let Some((name, value)) = parsed.filters().first() {
            if header.filters().contains_key(name) {
                return Err(invalid(format!(
                    "FILTER/{name} already exists in the VCF header"
                )));
            }
            header.filters_mut().insert(name.clone(), value.clone());
        } else if let Some((name, value)) = parsed.contigs().first() {
            if header.contigs().contains_key(name) {
                return Err(invalid(format!(
                    "contig/{name} already exists in the VCF header"
                )));
            }
            header.contigs_mut().insert(name.clone(), value.clone());
        }
    }
    Ok(())
}

fn read_chromosome_renames(path: Option<&Path>) -> Result<HashMap<String, String>> {
    let Some(path) = path else {
        return Ok(HashMap::new());
    };
    let mut mappings = HashMap::new();
    let mut targets = HashSet::new();
    for (old, new, line) in read_rename_pairs(path)? {
        if !valid_contig(&old) || !valid_contig(&new) {
            return Err(invalid(format!("invalid chromosome rename on line {line}")));
        }
        if mappings.insert(old, new.clone()).is_some() {
            return Err(invalid(format!(
                "chromosome rename map repeats a source on line {line}"
            )));
        }
        if !targets.insert(new.clone()) {
            return Err(invalid(format!(
                "chromosome rename map repeats target {new:?}"
            )));
        }
    }
    Ok(mappings)
}

fn read_annotation_renames(path: Option<&Path>) -> Result<Renames> {
    let Some(path) = path else {
        return Ok(Renames::default());
    };
    let mut renames = Renames::default();
    let mut targets = HashSet::new();
    for (old, new, line) in read_rename_pairs(path)? {
        let (namespace, old) = typed_name(&old)
            .map_err(|_| invalid(format!("invalid annotation rename source on line {line}")))?;
        let new = match new.split_once('/') {
            Some(_) => {
                let (new_namespace, new) = typed_name(&new).map_err(|_| {
                    invalid(format!("invalid annotation rename target on line {line}"))
                })?;
                if namespace != new_namespace {
                    return Err(invalid(format!(
                        "annotation rename crosses namespaces on line {line}"
                    )));
                }
                new
            }
            None if valid_tag(&new) => new,
            None => {
                return Err(invalid(format!(
                    "invalid annotation rename target on line {line}"
                )));
            }
        };
        if !targets.insert((namespace, new.clone())) {
            return Err(invalid(format!("multiple annotation renames target {new}")));
        }
        let map = match namespace {
            Namespace::Info => &mut renames.infos,
            Namespace::Format => &mut renames.formats,
            Namespace::Filter => &mut renames.filters,
        };
        if map.insert(old, new).is_some() {
            return Err(invalid(format!(
                "annotation rename map repeats a source on line {line}"
            )));
        }
    }
    Ok(renames)
}

fn read_rename_pairs(path: &Path) -> Result<Vec<(String, String, usize)>> {
    let content = fs::read_to_string(path)
        .rs_with_context(|| format!("reading rename map {}", path.display()))?;
    let mut mappings = Vec::new();
    for (index, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() != 2 {
            return Err(invalid(format!(
                "rename map {} line {} must contain exactly two fields",
                path.display(),
                index + 1
            )));
        }
        mappings.push((fields[0].to_owned(), fields[1].to_owned(), index + 1));
    }
    if mappings.is_empty() {
        return Err(invalid(format!("rename map is empty: {}", path.display())));
    }
    Ok(mappings)
}

fn validate_removals(header: &Header, removals: &[Removal]) -> Result<()> {
    for removal in removals {
        match removal {
            Removal::Filters(selection) => validate_selection(selection, "FILTER", |name| {
                header.filters().contains_key(name)
            })?,
            Removal::Info(selection) => {
                validate_selection(selection, "INFO", |name| header.infos().contains_key(name))?
            }
            Removal::Format(selection) => validate_selection(selection, "FORMAT", |name| {
                header.formats().contains_key(name)
            })?,
            Removal::Id | Removal::Qual => {}
        }
    }
    Ok(())
}

fn validate_selection(
    selection: &Selection<String>,
    kind: &str,
    contains: impl Fn(&str) -> bool,
) -> Result<()> {
    for value in &selection.values {
        if !contains(value) {
            return Err(invalid(format!(
                "{kind}/{value} is not defined in the VCF header"
            )));
        }
    }
    Ok(())
}

fn remove_header_definitions(header: &mut Header, removals: &[Removal]) {
    for removal in removals {
        match removal {
            Removal::Filters(selection) => header
                .filters_mut()
                .retain(|name, _| !selection.removes(name)),
            Removal::Info(selection) => header
                .infos_mut()
                .retain(|name, _| !selection.removes(name)),
            Removal::Format(selection) => header.formats_mut().retain(|name, _| {
                !selection.removes(name)
                    || (!selection.keep && selection.values.is_empty() && name == "GT")
            }),
            Removal::Id | Removal::Qual => {}
        }
    }
}

fn validate_renames(header: &Header, renames: &Renames) -> Result<()> {
    validate_map(&renames.chromosomes, "contig", |name| {
        header.contigs().contains_key(name)
    })?;
    validate_map(&renames.infos, "INFO", |name| {
        header.infos().contains_key(name)
    })?;
    validate_map(&renames.formats, "FORMAT", |name| {
        header.formats().contains_key(name)
    })?;
    validate_map(&renames.filters, "FILTER", |name| {
        header.filters().contains_key(name)
    })
}

fn validate_map(
    renames: &HashMap<String, String>,
    kind: &str,
    contains: impl Fn(&str) -> bool,
) -> Result<()> {
    for (old, new) in renames {
        if old == new {
            return Err(invalid(format!("{kind} rename does not change {old}")));
        }
        if !contains(old) {
            return Err(invalid(format!(
                "{kind}/{old} is not defined in the VCF header"
            )));
        }
        if contains(new) {
            return Err(invalid(format!(
                "{kind}/{new} already exists in the VCF header"
            )));
        }
    }
    Ok(())
}

fn rename_header_maps(header: &mut Header, renames: &Renames) {
    macro_rules! rename {
        ($map:expr, $renames:expr) => {
            for (old, new) in $renames {
                let index = $map.get_index_of(old).expect("validated rename source");
                let (_, value) = $map
                    .shift_remove_entry(old)
                    .expect("validated rename source");
                let replaced = $map.shift_insert(index, new.clone(), value);
                debug_assert!(replaced.is_none());
            }
        };
    }
    rename!(header.contigs_mut(), &renames.chromosomes);
    rename!(header.infos_mut(), &renames.infos);
    rename!(header.formats_mut(), &renames.formats);
    rename!(header.filters_mut(), &renames.filters);
}

fn apply_removals(record: &mut RecordBuf, removals: &[Removal]) -> bool {
    let mut changed = false;
    for removal in removals {
        match removal {
            Removal::Id => {
                changed |= !record.ids().as_ref().is_empty();
                record.ids_mut().as_mut().clear();
            }
            Removal::Qual => changed |= record.quality_score_mut().take().is_some(),
            Removal::Filters(selection) => {
                let filters = record.filters_mut().as_mut();
                let before = filters.len();
                filters.retain(|name| !selection.removes(name));
                changed |= filters.len() != before;
            }
            Removal::Info(selection) => {
                let info = record.info_mut().as_mut();
                let before = info.len();
                info.retain(|name, _| !selection.removes(name));
                changed |= info.len() != before;
            }
            Removal::Format(selection) => {
                let keys = record.samples().keys().as_ref();
                let retained: Vec<_> = keys
                    .iter()
                    .enumerate()
                    .filter(|(_, name)| {
                        !selection.removes(name)
                            || (!selection.keep
                                && selection.values.is_empty()
                                && name.as_str() == "GT")
                    })
                    .map(|(index, name)| (index, name.clone()))
                    .collect();
                if retained.len() != keys.len() {
                    let keys: Keys = retained.iter().map(|(_, name)| name.clone()).collect();
                    let values = record
                        .samples()
                        .values()
                        .map(|sample| {
                            retained
                                .iter()
                                .map(|(index, _)| {
                                    sample.values().get(*index).cloned().unwrap_or(None)
                                })
                                .collect()
                        })
                        .collect();
                    *record.samples_mut() = Samples::new(keys, values);
                    changed = true;
                }
            }
        }
    }
    changed
}

fn apply_renames(record: &mut RecordBuf, renames: &Renames) -> Result<bool> {
    let mut changed = false;
    if let Some(new) = renames.chromosomes.get(record.reference_sequence_name()) {
        *record.reference_sequence_name_mut() = new.clone();
        changed = true;
    }

    macro_rules! rename_map {
        ($map:expr, $renames:expr) => {
            for (old, new) in $renames {
                if let Some(index) = $map.get_index_of(old) {
                    if $map.contains_key(new) {
                        return Err(invalid_record(format!(
                            "annotation rename target {new} is already present"
                        )));
                    }
                    let (_, value) = $map.shift_remove_entry(old).expect("present record key");
                    let replaced = $map.shift_insert(index, new.clone(), value);
                    debug_assert!(replaced.is_none());
                    changed = true;
                }
            }
        };
    }
    rename_map!(record.info_mut().as_mut(), &renames.infos);

    for (old, new) in &renames.formats {
        let keys = record.samples_mut().keys_mut().as_mut();
        if let Some(index) = keys.get_index_of(old) {
            if keys.contains(new) {
                return Err(invalid_record(format!(
                    "FORMAT rename target {new} is already present"
                )));
            }
            keys.shift_remove(old);
            let inserted = keys.shift_insert(index, new.clone());
            debug_assert!(inserted);
            changed = true;
        }
    }

    for (old, new) in &renames.filters {
        let filters = record.filters_mut().as_mut();
        if let Some(index) = filters.get_index_of(old) {
            if filters.contains(new) {
                return Err(invalid_record(format!(
                    "FILTER rename target {new} is already present"
                )));
            }
            filters.shift_remove(old);
            let inserted = filters.shift_insert(index, new.clone());
            debug_assert!(inserted);
            changed = true;
        }
    }
    Ok(changed)
}

fn valid_tag(tag: &str) -> bool {
    !tag.is_empty()
        && tag
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.' | b'-'))
}

fn valid_contig(name: &str) -> bool {
    !name.is_empty()
        && name.bytes().all(|byte| {
            !byte.is_ascii_whitespace() && !matches!(byte, b',' | b'<' | b'>' | b'[' | b']')
        })
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::ConfigError(message.into())
}

fn invalid_record(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use std::fs;

    use noodles_vcf::{self as vcf, variant::RecordBuf};

    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.3\n\
##INFO=<ID=OLD,Number=1,Type=Integer,Description=\"old\">\n\
##INFO=<ID=KEEP,Number=1,Type=String,Description=\"keep\">\n\
##FILTER=<ID=q10,Description=\"q10\">\n\
##FILTER=<ID=s50,Description=\"s50\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"x\">\n\
##FORMAT=<ID=Y,Number=1,Type=Integer,Description=\"y\">\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";

    fn fixture() -> (vcf::Header, RecordBuf) {
        let header: vcf::Header = HEADER.parse().unwrap();
        let mut reader = vcf::io::Reader::new(
            b"chr1\t10\trs1\tA\tC\t12\tq10;s50\tOLD=1;KEEP=x\tGT:X:Y\t0/1:4:7\t1/1:5:8\n"
                .as_slice(),
        );
        let raw = reader.records().next().unwrap().unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        (header, record)
    }

    fn options(remove: Option<&str>) -> HeaderOptions {
        HeaderOptions {
            remove: remove.map(str::to_owned),
            ..HeaderOptions::default()
        }
    }

    fn bind_renames(
        header: &vcf::Header,
        chromosomes: &str,
        annotations: &str,
    ) -> rsomics_common::Result<HeaderPlan> {
        let directory = tempfile::tempdir().unwrap();
        let chromosome_path = directory.path().join("chromosomes.tsv");
        let annotation_path = directory.path().join("annotations.tsv");
        fs::write(&chromosome_path, chromosomes).unwrap();
        fs::write(&annotation_path, annotations).unwrap();
        HeaderPlan::bind(
            header,
            HeaderOptions {
                rename_chromosomes: Some(chromosome_path),
                rename_annotations: Some(annotation_path),
                ..HeaderOptions::default()
            },
        )
    }

    #[test]
    fn removes_fixed_fields_and_whole_namespaces() {
        let (mut header, mut record) = fixture();
        let plan = HeaderPlan::bind(&header, options(Some("ID,QUAL,FILTER,INFO,FORMAT"))).unwrap();

        plan.prepare(&mut header).unwrap();
        assert!(plan.apply(&mut record).unwrap());

        assert!(record.ids().as_ref().is_empty());
        assert!(record.quality_score().is_none());
        assert!(record.filters().as_ref().is_empty());
        assert!(record.info().as_ref().is_empty());
        assert_eq!(record.samples().keys().as_ref().len(), 1);
        assert!(record.samples().keys().as_ref().contains("GT"));
        assert!(header.filters().is_empty());
        assert!(header.infos().is_empty());
        assert_eq!(header.formats().len(), 1);
        assert!(header.formats().contains_key("GT"));
    }

    #[test]
    fn removes_selected_tags_and_retains_complements() {
        let (mut header, mut record) = fixture();
        let plan = HeaderPlan::bind(
            &header,
            options(Some("INFO/OLD,^FORMAT/GT,FORMAT/Y,^FILTER/q10")),
        )
        .unwrap();

        plan.prepare(&mut header).unwrap();
        plan.apply(&mut record).unwrap();

        assert!(!record.info().as_ref().contains_key("OLD"));
        assert!(record.info().as_ref().contains_key("KEEP"));
        let keys = record.samples().keys().as_ref();
        assert_eq!(
            keys.iter().map(String::as_str).collect::<Vec<_>>(),
            ["GT", "Y"]
        );
        assert_eq!(
            record
                .filters()
                .as_ref()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["q10"]
        );
        assert_eq!(
            header
                .formats()
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["GT", "Y"]
        );
        assert_eq!(
            header
                .filters()
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["q10"]
        );
    }

    #[test]
    fn removes_gt_when_requested_explicitly() {
        let (mut header, mut record) = fixture();
        let plan = HeaderPlan::bind(&header, options(Some("FORMAT/GT"))).unwrap();

        plan.prepare(&mut header).unwrap();
        plan.apply(&mut record).unwrap();

        assert!(!header.formats().contains_key("GT"));
        assert!(!record.samples().keys().as_ref().contains("GT"));
    }

    #[test]
    fn appends_supported_header_definitions() {
        let (mut header, _) = fixture();
        let plan = HeaderPlan::bind(
            &header,
            HeaderOptions {
                appended: vec![
                    "##INFO=<ID=NEW,Number=1,Type=Float,Description=\"new\">".into(),
                    "##FORMAT=<ID=Z,Number=1,Type=String,Description=\"z\">".into(),
                    "##FILTER=<ID=LowQ,Description=\"low\">".into(),
                    "##contig=<ID=chr2,length=200>".into(),
                ],
                ..HeaderOptions::default()
            },
        )
        .unwrap();

        plan.prepare(&mut header).unwrap();

        assert!(header.infos().contains_key("NEW"));
        assert!(header.formats().contains_key("Z"));
        assert!(header.filters().contains_key("LowQ"));
        assert!(header.contigs().contains_key("chr2"));
    }

    #[test]
    fn renames_header_and_record_keys_together() {
        let (mut header, mut record) = fixture();
        let plan = bind_renames(
            &header,
            "chr1\t1\n",
            "INFO/OLD\tNEW\nFORMAT/X\tNEW\nFILTER/q10\tLowQ\n",
        )
        .unwrap();

        plan.prepare(&mut header).unwrap();
        plan.apply(&mut record).unwrap();

        assert!(header.contigs().contains_key("1"));
        assert!(header.infos().contains_key("NEW"));
        assert!(header.formats().contains_key("NEW"));
        assert!(header.filters().contains_key("LowQ"));
        assert_eq!(record.reference_sequence_name(), "1");
        assert!(record.info().as_ref().contains_key("NEW"));
        assert!(record.samples().keys().as_ref().contains("NEW"));
        assert!(record.filters().as_ref().contains("LowQ"));
    }

    #[test]
    fn rejects_invalid_appended_lines_and_missing_names() {
        let (header, _) = fixture();

        for line in [
            "##source=bad",
            "##INFO=<ID=OLD,Number=1,Type=Integer,Description=\"duplicate\">",
        ] {
            let result = HeaderPlan::bind(
                &header,
                HeaderOptions {
                    appended: vec![line.into()],
                    ..HeaderOptions::default()
                },
            );
            assert!(result.is_err(), "{line}");
        }

        assert!(HeaderPlan::bind(&header, options(Some("INFO/MISSING"))).is_err());
        assert!(bind_renames(&header, "missing\t1\n", "INFO/OLD\tNEW\n").is_err());
    }

    #[test]
    fn rejects_rename_collisions_and_duplicate_targets() {
        let (header, _) = fixture();

        for annotations in [
            "INFO/OLD\tKEEP\n",
            "INFO/OLD\tNEW\nINFO/KEEP\tNEW\n",
            "FORMAT/X\tINFO/Z\n",
        ] {
            assert!(bind_renames(&header, "chr1\t1\n", annotations).is_err());
        }

        assert!(bind_renames(&header, "chr1\t1\nchr1\t2\n", "INFO/OLD\tNEW\n").is_err());
    }

    #[test]
    fn rejects_record_level_rename_collisions() {
        let (header, mut record) = fixture();
        let plan = bind_renames(&header, "chr1\t1\n", "INFO/OLD\tNEW\n").unwrap();
        let value = record.info().as_ref().get("OLD").cloned().unwrap();
        record.info_mut().as_mut().insert("NEW".into(), value);

        assert!(plan.apply(&mut record).is_err());
    }
}
