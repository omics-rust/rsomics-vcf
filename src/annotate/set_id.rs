use noodles_vcf::{Header, variant::RecordBuf};
use rsomics_common::{Result, RsomicsError};

use crate::{format::HeaderTypes, query_format::SiteFormat};

pub(crate) struct IdPlan {
    only_missing: bool,
    format: SiteFormat,
    scratch: Vec<u8>,
}

impl IdPlan {
    pub(crate) fn bind(source: Option<&str>, schema: &HeaderTypes) -> Result<Option<Self>> {
        let Some(source) = source else {
            return Ok(None);
        };
        let (only_missing, source) = source
            .strip_prefix('+')
            .map_or((false, source), |source| (true, source));
        if source.is_empty() {
            return Err(invalid("set-id format must not be empty"));
        }
        Ok(Some(Self {
            only_missing,
            format: SiteFormat::bind(source, schema)?,
            scratch: Vec::new(),
        }))
    }

    pub(crate) fn apply(&mut self, header: &Header, record: &mut RecordBuf) -> Result<bool> {
        if self.only_missing && !record.ids().as_ref().is_empty() {
            return Ok(false);
        }
        self.format.render(header, record, &mut self.scratch)?;
        validate_rendered(&self.scratch)?;

        let ids = record.ids_mut().as_mut();
        let previous: Vec<_> = ids.iter().cloned().collect();
        ids.clear();
        if self.scratch != b"." {
            ids.extend(
                self.scratch
                    .split(|byte| *byte == b';')
                    .map(|value| String::from_utf8(value.to_vec()).expect("validated UTF-8 ID")),
            );
        }
        Ok(ids.iter().ne(previous.iter()))
    }
}

fn validate_rendered(value: &[u8]) -> Result<()> {
    if value.is_empty() {
        return Err(invalid("set-id format rendered an empty ID"));
    }
    let text = std::str::from_utf8(value)
        .map_err(|_| invalid("set-id format rendered an ID that is not UTF-8"))?;
    if text == "." {
        return Ok(());
    }
    if text
        .split(';')
        .any(|id| id.is_empty() || id == "." || id.bytes().any(|byte| byte.is_ascii_whitespace()))
    {
        return Err(invalid(format!(
            "set-id format rendered an invalid ID: {text:?}"
        )));
    }
    Ok(())
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::ConfigError(message.into())
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{self as vcf, variant::RecordBuf};

    use crate::format::HeaderTypes;

    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.3\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
##INFO=<ID=S,Number=1,Type=String,Description=\"string\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn record(id: &str) -> (vcf::Header, HeaderTypes, RecordBuf) {
        let header: vcf::Header = HEADER.parse().unwrap();
        let schema = HeaderTypes::parse(HEADER.as_bytes()).unwrap();
        let line = format!("chr1\t10\t{id}\tA\tC,G\t12\tPASS\tDP=7;S=bad%20value\n");
        let mut reader = vcf::io::Reader::new(line.as_bytes());
        let raw = reader.records().next().unwrap().unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        (header, schema, record)
    }

    #[test]
    fn replaces_ids_from_site_fields() {
        let (header, schema, mut record) = record("old");
        let mut plan = IdPlan::bind(Some(r"%CHROM\_%POS\_%REF\_%FIRST_ALT"), &schema)
            .unwrap()
            .unwrap();

        assert!(plan.apply(&header, &mut record).unwrap());
        assert_eq!(
            record
                .ids()
                .as_ref()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["chr1_10_A_C"]
        );
    }

    #[test]
    fn plus_sets_only_missing_ids() {
        let (header, schema, mut present) = record("old");
        let (_, _, mut missing) = record(".");
        let mut plan = IdPlan::bind(Some(r"+%CHROM\_%POS"), &schema)
            .unwrap()
            .unwrap();

        assert!(!plan.apply(&header, &mut present).unwrap());
        assert!(plan.apply(&header, &mut missing).unwrap());
        assert!(present.ids().as_ref().contains("old"));
        assert!(missing.ids().as_ref().contains("chr1_10"));
    }

    #[test]
    fn rejects_invalid_rendered_ids() {
        let (header, schema, mut record) = record("old");
        let plan = IdPlan::bind(Some("%INFO/MISSING"), &schema);
        assert!(plan.is_err());

        record.info_mut().as_mut().insert(
            "S".into(),
            Some(noodles_vcf::variant::record_buf::info::field::Value::String("bad value".into())),
        );
        let mut plan = IdPlan::bind(Some("%INFO/S"), &schema).unwrap().unwrap();
        assert!(plan.apply(&header, &mut record).is_err());
    }
}
