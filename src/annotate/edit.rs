mod samples;

use std::collections::HashSet;

use noodles_vcf::{
    Header,
    header::record::value::map::info::{Number, Type},
    variant::{
        RecordBuf,
        record_buf::{
            Filters, Ids,
            info::field::{Value, value::Array},
        },
    },
};
use rsomics_common::{Result, RsomicsError};

pub(crate) use self::samples::SampleSelection;
use self::samples::{FormatPlan, prepare_definition as prepare_format_definition};

use super::{
    columns::{BoundColumns, Destination, SourceField, SourceKind, Transfer, WriteMode},
    matching::Matched,
    source::Payload,
};
use crate::norm::cardinality::{combinations, genotype_index, infer_ploidy, visit_genotypes};

#[derive(Clone)]
struct InfoPlan {
    number: Number,
    ty: Type,
    context: String,
}

#[derive(Clone)]
struct BoundTransfer {
    transfer: Transfer,
    info: Option<InfoPlan>,
    format: Option<FormatPlan>,
}

pub(crate) struct Editor {
    source_kind: SourceKind,
    transfers: Vec<BoundTransfer>,
    samples: Option<SampleSelection>,
}

#[derive(Clone, Debug, PartialEq)]
enum State<T> {
    Absent,
    Missing,
    Remove,
    Value(T),
}

impl Editor {
    pub(crate) fn bind(
        target: &mut Header,
        source: Option<&Header>,
        columns: BoundColumns,
    ) -> Result<Self> {
        let source_kind = columns.source_kind();
        let mut prepared = target.clone();
        let mut transfers = Vec::new();
        let mut destinations = HashSet::new();
        for transfer in columns.spec().transfers() {
            match (&transfer.source, &transfer.destination) {
                (SourceField::AllInfo, Destination::AllInfo) => {
                    let source =
                        source.ok_or_else(|| invalid("INFO transfer requires a source header"))?;
                    for key in source.infos().keys() {
                        if !destinations.insert(format!("INFO/{key}")) {
                            return Err(invalid(format!(
                                "annotation destination is assigned more than once: INFO/{key}"
                            )));
                        }
                        let transfer = Transfer {
                            source: SourceField::Info(key.clone()),
                            destination: Destination::Info(key.clone()),
                            mode: transfer.mode,
                        };
                        transfers.push(bind_transfer(
                            &mut prepared,
                            Some(source),
                            transfer,
                            source_kind,
                        )?);
                    }
                }
                (SourceField::AllFormat, Destination::AllFormat) => {
                    let source = source
                        .ok_or_else(|| invalid("FORMAT transfer requires a source header"))?;
                    for key in source.formats().keys().filter(|key| *key != "GT") {
                        if !destinations.insert(format!("FORMAT/{key}")) {
                            return Err(invalid(format!(
                                "annotation destination is assigned more than once: FORMAT/{key}"
                            )));
                        }
                        let transfer = Transfer {
                            source: SourceField::Format(key.clone()),
                            destination: Destination::Format(key.clone()),
                            mode: transfer.mode,
                        };
                        transfers.push(bind_transfer(
                            &mut prepared,
                            Some(source),
                            transfer,
                            source_kind,
                        )?);
                    }
                }
                _ => {
                    let destination = destination_key(&transfer.destination);
                    if !destinations.insert(destination.clone()) {
                        return Err(invalid(format!(
                            "annotation destination is assigned more than once: {destination}"
                        )));
                    }
                    transfers.push(bind_transfer(
                        &mut prepared,
                        source,
                        transfer.clone(),
                        source_kind,
                    )?);
                }
            }
        }
        let samples = if transfers.iter().any(|bound| bound.format.is_some()) {
            let source =
                source.ok_or_else(|| invalid("FORMAT transfer requires a source header"))?;
            Some(SampleSelection::bind(source, &prepared, None, None)?)
        } else {
            None
        };
        *target = prepared;
        Ok(Self {
            source_kind,
            transfers,
            samples,
        })
    }

    pub(crate) fn apply_info(
        &self,
        header: &Header,
        matched: &Matched<'_>,
        target: &mut RecordBuf,
    ) -> Result<bool> {
        let mut changed = false;
        for bound in &self.transfers {
            changed |= match &bound.transfer.destination {
                Destination::Id => self.apply_ids(header, bound, matched, target)?,
                Destination::Qual => self.apply_qual(bound, matched, target)?,
                Destination::Filter => self.apply_filters(header, bound, matched, target)?,
                Destination::Info(key) => self.apply_one_info(bound, key, matched, target)?,
                Destination::Format(_) | Destination::AllInfo | Destination::AllFormat => false,
            };
        }
        Ok(changed)
    }

    fn apply_ids(
        &self,
        _header: &Header,
        bound: &BoundTransfer,
        matched: &Matched<'_>,
        target: &mut RecordBuf,
    ) -> Result<bool> {
        let source = match source_bytes(&bound.transfer.source, matched)? {
            Some(Some(raw)) => State::Value(parse_ids(raw, matched.source.serial)?),
            Some(None) => State::Missing,
            None => match &matched.source.payload {
                Payload::Variant(record) => {
                    if record.ids().as_ref().is_empty() {
                        State::Missing
                    } else {
                        State::Value(record.ids().clone())
                    }
                }
                Payload::Tabular(_) => State::Absent,
            },
        };
        let current = if target.ids().as_ref().is_empty() {
            State::Missing
        } else {
            State::Value(target.ids().clone())
        };
        let next = choose(
            bound.transfer.mode,
            source,
            current.clone(),
            self.source_kind,
        )?;
        if next == current {
            return Ok(false);
        }
        *target.ids_mut() = match next {
            State::Value(ids) => ids,
            State::Absent | State::Missing | State::Remove => Ids::default(),
        };
        Ok(true)
    }

    fn apply_qual(
        &self,
        bound: &BoundTransfer,
        matched: &Matched<'_>,
        target: &mut RecordBuf,
    ) -> Result<bool> {
        let source = match source_bytes(&bound.transfer.source, matched)? {
            Some(Some(raw)) => State::Value(parse_quality(raw, matched.source.serial)?),
            Some(None) => State::Missing,
            None => match &matched.source.payload {
                Payload::Variant(record) => {
                    record.quality_score().map_or(State::Missing, State::Value)
                }
                Payload::Tabular(_) => State::Absent,
            },
        };
        let current = target.quality_score().map_or(State::Missing, State::Value);
        let next = choose(
            bound.transfer.mode,
            source,
            current.clone(),
            self.source_kind,
        )?;
        if next == current {
            return Ok(false);
        }
        *target.quality_score_mut() = match next {
            State::Value(value) => Some(value),
            State::Absent | State::Missing | State::Remove => None,
        };
        Ok(true)
    }

    fn apply_filters(
        &self,
        header: &Header,
        bound: &BoundTransfer,
        matched: &Matched<'_>,
        target: &mut RecordBuf,
    ) -> Result<bool> {
        let source = match source_bytes(&bound.transfer.source, matched)? {
            Some(Some(raw)) => State::Value(parse_filters(raw, header, matched.source.serial)?),
            Some(None) => State::Missing,
            None => match &matched.source.payload {
                Payload::Variant(record) => {
                    if record.filters().as_ref().is_empty() {
                        State::Missing
                    } else {
                        State::Value(record.filters().clone())
                    }
                }
                Payload::Tabular(_) => State::Absent,
            },
        };
        let current = if target.filters().as_ref().is_empty() {
            State::Missing
        } else {
            State::Value(target.filters().clone())
        };
        let next = choose(
            bound.transfer.mode,
            source,
            current.clone(),
            self.source_kind,
        )?;
        if next == current {
            return Ok(false);
        }
        *target.filters_mut() = match next {
            State::Value(filters) => filters,
            State::Absent | State::Missing | State::Remove => Filters::default(),
        };
        Ok(true)
    }

    fn apply_one_info(
        &self,
        bound: &BoundTransfer,
        key: &str,
        matched: &Matched<'_>,
        target: &mut RecordBuf,
    ) -> Result<bool> {
        let plan = bound.info.as_ref().expect("INFO transfer has a schema");
        let source = match source_bytes(&bound.transfer.source, matched)? {
            Some(Some(raw)) if plan.ty == Type::Flag => {
                parse_flag(raw, &plan.context, matched.source.serial)?
            }
            Some(Some(raw)) => State::Value(parse_info(raw, plan, matched.source.serial)?),
            Some(None) if plan.ty == Type::Flag => State::Absent,
            Some(None) => State::Missing,
            None => match (&matched.source.payload, &bound.transfer.source) {
                (Payload::Variant(record), SourceField::Info(source_key)) => {
                    match record.info().get(source_key) {
                        None => State::Absent,
                        Some(None) => State::Missing,
                        Some(Some(value)) => State::Value(value.clone()),
                    }
                }
                _ => State::Absent,
            },
        };
        let source = match source {
            State::Value(value) => {
                validate_info(&value, plan, matched.source.alternates.len() + 1)?;
                State::Value(remap_info(
                    plan.number,
                    value,
                    matched,
                    target.alternate_bases().as_ref().len() + 1,
                    &plan.context,
                )?)
            }
            state => state,
        };
        let current = match target.info().get(key) {
            None => State::Absent,
            Some(None) => State::Missing,
            Some(Some(value)) => State::Value(value.clone()),
        };
        let next = if matches!(
            bound.transfer.mode,
            WriteMode::Append | WriteMode::AppendMissing
        ) {
            append_state(
                bound.transfer.mode,
                source,
                current.clone(),
                plan.ty,
                &plan.context,
            )?
        } else {
            choose(
                bound.transfer.mode,
                source,
                current.clone(),
                self.source_kind,
            )?
        };
        if let State::Value(value) = &next {
            validate_info(value, plan, target.alternate_bases().as_ref().len() + 1)?;
        }
        if next == current {
            return Ok(false);
        }
        match next {
            State::Absent => {
                target.info_mut().as_mut().shift_remove(key);
            }
            State::Remove => unreachable!("remove is resolved by the write policy"),
            State::Missing => {
                target.info_mut().insert(key.to_owned(), None);
            }
            State::Value(value) => {
                target.info_mut().insert(key.to_owned(), Some(value));
            }
        }
        Ok(true)
    }
}

fn bind_transfer(
    target: &mut Header,
    source: Option<&Header>,
    transfer: Transfer,
    source_kind: SourceKind,
) -> Result<BoundTransfer> {
    let compatible = matches!(
        (&transfer.source, &transfer.destination),
        (SourceField::Id, Destination::Id)
            | (SourceField::Qual, Destination::Qual)
            | (SourceField::Filter, Destination::Filter)
            | (SourceField::Info(_), Destination::Info(_))
            | (SourceField::Format(_), Destination::Format(_))
            | (SourceField::AllInfo, Destination::AllInfo)
            | (SourceField::AllFormat, Destination::AllFormat)
            | (SourceField::Tabular(_), _)
    );
    if !compatible {
        return Err(invalid(
            "annotation source and destination fields are incompatible",
        ));
    }
    if matches!(transfer.mode, WriteMode::Append | WriteMode::AppendMissing)
        && !matches!(transfer.destination, Destination::Info(_))
    {
        return Err(invalid(
            "annotation append mode is not supported for this destination",
        ));
    }
    if matches!(transfer.destination, Destination::Filter)
        && let Some(source) = source
    {
        for (key, definition) in source.filters() {
            match target.filters().get(key) {
                Some(existing) if existing != definition => {
                    return Err(invalid(format!(
                        "FILTER/{key} has incompatible definitions"
                    )));
                }
                Some(_) => {}
                None => {
                    target.filters_mut().insert(key.clone(), definition.clone());
                }
            }
        }
    }
    let info = match (&transfer.source, &transfer.destination) {
        (SourceField::Info(source_key), Destination::Info(target_key)) => {
            let source = source.ok_or_else(|| invalid("INFO transfer requires a source header"))?;
            let definition = source
                .infos()
                .get(source_key)
                .ok_or_else(|| invalid(format!("source header has no INFO/{source_key}")))?
                .clone();
            prepare_info_definition(target, target_key, &definition)?;
            Some(InfoPlan {
                number: definition.number(),
                ty: definition.ty(),
                context: format!("INFO/{target_key}"),
            })
        }
        (SourceField::Tabular(_), Destination::Info(target_key)) => {
            let mut definition = target
                .infos()
                .get(target_key)
                .ok_or_else(|| invalid(format!("target header has no INFO/{target_key}")))?
                .clone();
            if matches!(transfer.mode, WriteMode::Append | WriteMode::AppendMissing) {
                if definition.ty() == Type::Flag {
                    return Err(invalid(format!(
                        "INFO/{target_key} flags cannot be appended"
                    )));
                }
                *definition.number_mut() = Number::Unknown;
                target
                    .infos_mut()
                    .insert(target_key.clone(), definition.clone());
            }
            Some(InfoPlan {
                number: definition.number(),
                ty: definition.ty(),
                context: format!("INFO/{target_key}"),
            })
        }
        (_, Destination::Info(key)) => {
            return Err(invalid(format!(
                "INFO/{key} has an incompatible source field"
            )));
        }
        _ => None,
    };
    let format = match (&transfer.source, &transfer.destination) {
        (SourceField::Format(source_key), Destination::Format(target_key)) => {
            if (source_key == "GT") != (target_key == "GT") {
                return Err(invalid("FORMAT/GT cannot be renamed"));
            }
            let source =
                source.ok_or_else(|| invalid("FORMAT transfer requires a source header"))?;
            let definition = source
                .formats()
                .get(source_key)
                .ok_or_else(|| invalid(format!("source header has no FORMAT/{source_key}")))?
                .clone();
            prepare_format_definition(target, target_key, &definition)?;
            Some(FormatPlan {
                number: definition.number(),
                ty: definition.ty(),
                genotype: target_key == "GT",
                context: format!("FORMAT/{target_key}"),
            })
        }
        (_, Destination::Format(key)) => {
            return Err(invalid(format!(
                "FORMAT/{key} has an incompatible source field"
            )));
        }
        _ => None,
    };
    if source_kind == SourceKind::Variant
        && matches!(transfer.mode, WriteMode::Append | WriteMode::AppendMissing)
    {
        return Err(invalid("append mode requires a tabular annotation source"));
    }
    Ok(BoundTransfer {
        transfer,
        info,
        format,
    })
}

fn destination_key(destination: &Destination) -> String {
    match destination {
        Destination::Id => "ID".to_owned(),
        Destination::Qual => "QUAL".to_owned(),
        Destination::Filter => "FILTER".to_owned(),
        Destination::Info(key) => format!("INFO/{key}"),
        Destination::Format(key) => format!("FORMAT/{key}"),
        Destination::AllInfo => "INFO".to_owned(),
        Destination::AllFormat => "FORMAT".to_owned(),
    }
}

fn prepare_info_definition(
    target: &mut Header,
    key: &str,
    source: &noodles_vcf::header::record::value::Map<noodles_vcf::header::record::value::map::Info>,
) -> Result<()> {
    match target.infos().get(key) {
        Some(existing) if existing.number() != source.number() || existing.ty() != source.ty() => {
            Err(invalid(format!("INFO/{key} has an incompatible schema")))
        }
        Some(_) => Ok(()),
        None => {
            target.infos_mut().insert(key.to_owned(), source.clone());
            Ok(())
        }
    }
}

fn source_bytes<'a>(
    source: &SourceField,
    matched: &'a Matched<'_>,
) -> Result<Option<Option<&'a [u8]>>> {
    let SourceField::Tabular(index) = source else {
        return Ok(None);
    };
    let Payload::Tabular(fields) = &matched.source.payload else {
        return Err(invalid(
            "tabular annotation column is bound to a variant record",
        ));
    };
    fields
        .get(*index)
        .map(|value| value.as_deref())
        .map(Some)
        .ok_or_else(|| {
            invalid(format!(
                "annotation record {} is missing column {}",
                matched.source.serial,
                index + 1
            ))
        })
}

fn choose<T: Clone>(
    mode: WriteMode,
    source: State<T>,
    target: State<T>,
    source_kind: SourceKind,
) -> Result<State<T>> {
    if matches!(source, State::Remove) {
        return Ok(State::Absent);
    }
    let concrete = matches!(source, State::Value(_));
    let present = !matches!(source, State::Absent)
        && (!matches!(source, State::Missing) || source_kind == SourceKind::Variant);
    let output = match mode {
        WriteMode::Replace if present => source,
        WriteMode::ReplaceMissing if !matches!(source, State::Absent) => source,
        WriteMode::Add if !matches!(target, State::Value(_)) && concrete => source,
        WriteMode::AddMissing
            if !matches!(target, State::Value(_)) && !matches!(source, State::Absent) =>
        {
            source
        }
        WriteMode::ReplaceExisting if !matches!(target, State::Absent) && present => source,
        WriteMode::Append | WriteMode::AppendMissing => {
            return Err(invalid("append mode requires an appendable field"));
        }
        _ => target,
    };
    Ok(output)
}

fn append_state(
    mode: WriteMode,
    source: State<Value>,
    target: State<Value>,
    ty: Type,
    context: &str,
) -> Result<State<Value>> {
    let include_missing = mode == WriteMode::AppendMissing;
    match (target, source) {
        (State::Absent, State::Value(value)) => Ok(State::Value(value)),
        (State::Missing, State::Value(value)) if include_missing => {
            prepend_missing(Some(value), ty, context).map(State::Value)
        }
        (State::Missing, State::Value(value)) => Ok(State::Value(value)),
        (State::Absent, State::Missing) if include_missing => Ok(State::Missing),
        (State::Absent, State::Remove) => Ok(State::Absent),
        (State::Absent, _) => Ok(State::Absent),
        (State::Missing, State::Missing) if include_missing => {
            prepend_missing(None, ty, context).map(State::Value)
        }
        (State::Missing, State::Remove) => Ok(State::Absent),
        (State::Missing, _) => Ok(State::Missing),
        (State::Value(target), State::Value(source)) => {
            append_values(target, Some(source), context).map(State::Value)
        }
        (State::Value(target), State::Missing) if include_missing => {
            append_values(target, None, context).map(State::Value)
        }
        (State::Value(_), State::Remove) => Ok(State::Absent),
        (target, _) => Ok(target),
    }
}

fn prepend_missing(source: Option<Value>, ty: Type, context: &str) -> Result<Value> {
    macro_rules! prepend {
        ($variant:ident, $source:expr) => {{
            let mut output = vec![None];
            match $source {
                Some(Value::$variant(value)) => output.push(Some(value)),
                Some(Value::Array(Array::$variant(values))) => output.extend(values),
                None => output.push(None),
                Some(_) => return Err(invalid(format!("{context} append type mismatch"))),
            }
            Ok(Value::Array(Array::$variant(output)))
        }};
    }
    match ty {
        Type::Integer => prepend!(Integer, source),
        Type::Float => prepend!(Float, source),
        Type::Character => prepend!(Character, source),
        Type::String => prepend!(String, source),
        Type::Flag => Err(invalid(format!("{context} flags cannot be appended"))),
    }
}

fn append_values(target: Value, source: Option<Value>, context: &str) -> Result<Value> {
    macro_rules! append {
        ($variant:ident, $target:expr, $source:expr) => {{
            let mut values = match $target {
                Value::$variant(value) => vec![Some(value)],
                Value::Array(Array::$variant(values)) => values,
                _ => return Err(invalid(format!("{context} append type mismatch"))),
            };
            match $source {
                None => values.push(None),
                Some(Value::$variant(value)) => values.push(Some(value)),
                Some(Value::Array(Array::$variant(source))) => values.extend(source),
                Some(_) => return Err(invalid(format!("{context} append type mismatch"))),
            }
            Ok(Value::Array(Array::$variant(values)))
        }};
    }
    match target {
        value @ Value::Integer(_) | value @ Value::Array(Array::Integer(_)) => {
            append!(Integer, value, source)
        }
        value @ Value::Float(_) | value @ Value::Array(Array::Float(_)) => {
            append!(Float, value, source)
        }
        value @ Value::Character(_) | value @ Value::Array(Array::Character(_)) => {
            append!(Character, value, source)
        }
        value @ Value::String(_) | value @ Value::Array(Array::String(_)) => {
            append!(String, value, source)
        }
        Value::Flag => Err(invalid(format!("{context} flags cannot be appended"))),
    }
}

fn parse_ids(raw: &[u8], serial: u64) -> Result<Ids> {
    let raw = text(raw, serial, "ID")?;
    let values: Vec<_> = raw.split(';').collect();
    if values.iter().any(|value| value.is_empty()) {
        return Err(invalid(format!(
            "annotation record {serial} has an invalid ID"
        )));
    }
    Ok(values.into_iter().map(str::to_owned).collect())
}

fn parse_quality(raw: &[u8], serial: u64) -> Result<f32> {
    let value: f32 = text(raw, serial, "QUAL")?
        .parse()
        .map_err(|_| invalid(format!("annotation record {serial} has an invalid QUAL")))?;
    if !value.is_finite() || value < 0.0 {
        return Err(invalid(format!(
            "annotation record {serial} has an invalid QUAL"
        )));
    }
    Ok(value)
}

fn parse_filters(raw: &[u8], header: &Header, serial: u64) -> Result<Filters> {
    let raw = text(raw, serial, "FILTER")?;
    if raw == "PASS" {
        return Ok(Filters::pass());
    }
    let mut filters = Filters::default();
    for filter in raw.split(';') {
        if filter.is_empty() || !header.filters().contains_key(filter) {
            return Err(invalid(format!(
                "annotation record {serial} uses undefined FILTER/{filter}"
            )));
        }
        filters.as_mut().insert(filter.to_owned());
    }
    Ok(filters)
}

fn parse_info(raw: &[u8], plan: &InfoPlan, serial: u64) -> Result<Value> {
    let raw = text(raw, serial, &plan.context)?;
    if plan.ty == Type::Flag {
        return Err(invalid(format!(
            "annotation record {serial} has an invalid {} flag",
            plan.context
        )));
    }
    let array = !matches!(plan.number, Number::Count(1));
    macro_rules! parse_values {
        ($ty:ty, $scalar:ident, $array_variant:ident) => {{
            let values = raw
                .split(',')
                .map(|value| {
                    if value == "." {
                        Ok(None)
                    } else {
                        value.parse::<$ty>().map(Some).map_err(|_| {
                            invalid(format!(
                                "annotation record {serial} has an invalid {} value",
                                plan.context
                            ))
                        })
                    }
                })
                .collect::<Result<Vec<_>>>()?;
            if array {
                Ok(Value::Array(Array::$array_variant(values)))
            } else if values.len() == 1 {
                values[0].map(Value::$scalar).ok_or_else(|| {
                    invalid(format!(
                        "annotation record {serial} has a missing {} value",
                        plan.context
                    ))
                })
            } else {
                Err(invalid(format!(
                    "annotation record {serial} has too many {} values",
                    plan.context
                )))
            }
        }};
    }
    match plan.ty {
        Type::Integer => parse_values!(i32, Integer, Integer),
        Type::Float => parse_values!(f32, Float, Float),
        Type::Character => {
            let values = raw
                .split(',')
                .map(|value| {
                    if value == "." {
                        Ok(None)
                    } else {
                        let mut chars = value.chars();
                        match (chars.next(), chars.next()) {
                            (Some(value), None) => Ok(Some(value)),
                            _ => Err(invalid(format!(
                                "annotation record {serial} has an invalid {} value",
                                plan.context
                            ))),
                        }
                    }
                })
                .collect::<Result<Vec<_>>>()?;
            if array {
                Ok(Value::Array(Array::Character(values)))
            } else if values.len() == 1 {
                values[0].map(Value::Character).ok_or_else(|| {
                    invalid(format!(
                        "annotation record {serial} has a missing {} value",
                        plan.context
                    ))
                })
            } else {
                Err(invalid(format!(
                    "annotation record {serial} has too many {} values",
                    plan.context
                )))
            }
        }
        Type::String => {
            let values = raw
                .split(',')
                .map(|value| (value != ".").then(|| value.to_owned()))
                .collect::<Vec<_>>();
            if array {
                Ok(Value::Array(Array::String(values)))
            } else if values.len() == 1 {
                values
                    .into_iter()
                    .next()
                    .flatten()
                    .map(Value::String)
                    .ok_or_else(|| {
                        invalid(format!(
                            "annotation record {serial} has a missing {} value",
                            plan.context
                        ))
                    })
            } else {
                Err(invalid(format!(
                    "annotation record {serial} has too many {} values",
                    plan.context
                )))
            }
        }
        Type::Flag => unreachable!(),
    }
}

fn parse_flag(raw: &[u8], context: &str, serial: u64) -> Result<State<Value>> {
    match text(raw, serial, context)? {
        "0" => Ok(State::Remove),
        "1" => Ok(State::Value(Value::Flag)),
        _ => Err(invalid(format!(
            "annotation record {serial} has an invalid {context} flag"
        ))),
    }
}

fn validate_info(value: &Value, plan: &InfoPlan, alleles: usize) -> Result<()> {
    let (ty, count) = match value {
        Value::Integer(_) => (Type::Integer, 1),
        Value::Float(_) => (Type::Float, 1),
        Value::Flag => (Type::Flag, 0),
        Value::Character(_) => (Type::Character, 1),
        Value::String(_) => (Type::String, 1),
        Value::Array(Array::Integer(values)) => (Type::Integer, values.len()),
        Value::Array(Array::Float(values)) => (Type::Float, values.len()),
        Value::Array(Array::Character(values)) => (Type::Character, values.len()),
        Value::Array(Array::String(values)) => (Type::String, values.len()),
    };
    if ty != plan.ty {
        return Err(invalid(format!(
            "{} has type {ty}, expected {}",
            plan.context, plan.ty
        )));
    }
    let valid = match plan.number {
        Number::Count(expected) => count == expected,
        Number::AlternateBases => count == alleles - 1,
        Number::ReferenceAlternateBases => count == alleles,
        Number::Samples => infer_ploidy(alleles, count).is_some(),
        Number::Unknown => true,
    };
    if !valid {
        return Err(invalid(format!(
            "{} has invalid cardinality {count}",
            plan.context
        )));
    }
    Ok(())
}

fn remap_info(
    number: Number,
    value: Value,
    matched: &Matched<'_>,
    target_alleles: usize,
    context: &str,
) -> Result<Value> {
    if !matches!(number, Number::A | Number::R | Number::G) {
        return Ok(value);
    }
    macro_rules! remap {
        ($variant:ident, $values:expr) => {
            remap_array(number, $values, matched, target_alleles, context)
                .map(Array::$variant)
                .map(Value::Array)
        };
    }
    match value {
        Value::Array(Array::Integer(values)) => remap!(Integer, values),
        Value::Array(Array::Float(values)) => remap!(Float, values),
        Value::Array(Array::Character(values)) => remap!(Character, values),
        Value::Array(Array::String(values)) => remap!(String, values),
        _ => Err(invalid(format!(
            "{context} allele-indexed value must be an array"
        ))),
    }
}

fn remap_array<T: Clone + PartialEq>(
    number: Number,
    source: Vec<Option<T>>,
    matched: &Matched<'_>,
    target_alleles: usize,
    context: &str,
) -> Result<Vec<Option<T>>> {
    let mapping = &matched.allele_map;
    if mapping.is_empty() || mapping[0] != Some(0) {
        return Err(invalid(format!(
            "{context} has no reference allele mapping"
        )));
    }
    if mapping
        .iter()
        .flatten()
        .any(|&index| index >= target_alleles)
    {
        return Err(invalid(format!(
            "{context} allele mapping exceeds the target"
        )));
    }
    match number {
        Number::AlternateBases | Number::ReferenceAlternateBases => {
            let offset = usize::from(number == Number::AlternateBases);
            let expected = mapping.len() - offset;
            if source.len() != expected {
                return Err(invalid(format!(
                    "{context} has {} values, expected {expected}",
                    source.len()
                )));
            }
            let target_count = target_alleles - offset;
            let mut output = vec![None; target_count];
            let mut assigned = vec![false; target_count];
            for (source_allele, &target_allele) in mapping.iter().enumerate().skip(offset) {
                let Some(target_allele) = target_allele else {
                    continue;
                };
                let source_index = source_allele - offset;
                let target_index = target_allele.checked_sub(offset).ok_or_else(|| {
                    invalid(format!("{context} alternate allele maps to the reference"))
                })?;
                assign(
                    &mut output,
                    &mut assigned,
                    target_index,
                    source[source_index].clone(),
                    context,
                )?;
            }
            Ok(output)
        }
        Number::Samples => {
            let ploidy = infer_ploidy(mapping.len(), source.len())
                .ok_or_else(|| invalid(format!("{context} genotype cardinality is ambiguous")))?;
            let genotype_width = target_alleles
                .checked_add(ploidy - 1)
                .ok_or_else(|| invalid(format!("{context} genotype cardinality overflows")))?;
            let target_count = combinations(genotype_width, ploidy)
                .ok_or_else(|| invalid(format!("{context} genotype cardinality overflows")))?;
            let mut output = vec![None; target_count];
            let mut assigned = vec![false; target_count];
            let mut error = None;
            visit_genotypes(mapping.len(), ploidy, |genotype| {
                if error.is_some() {
                    return;
                }
                let Some(mut target) = genotype
                    .iter()
                    .map(|&allele| mapping[allele])
                    .collect::<Option<Vec<_>>>()
                else {
                    return;
                };
                target.sort_unstable();
                let source_index = genotype_index(genotype).expect("enumerated genotype is valid");
                let target_index = genotype_index(&target).expect("mapped genotype is valid");
                if let Err(current) = assign(
                    &mut output,
                    &mut assigned,
                    target_index,
                    source[source_index].clone(),
                    context,
                ) {
                    error = Some(current);
                }
            });
            if let Some(error) = error {
                Err(error)
            } else {
                Ok(output)
            }
        }
        Number::Count(_) | Number::Unknown => unreachable!(),
    }
}

fn assign<T: PartialEq>(
    output: &mut [Option<T>],
    assigned: &mut [bool],
    index: usize,
    value: Option<T>,
    context: &str,
) -> Result<()> {
    if assigned[index] && output[index] != value {
        return Err(invalid(format!("{context} allele mapping is ambiguous")));
    }
    output[index] = value;
    assigned[index] = true;
    Ok(())
}

fn text<'a>(raw: &'a [u8], serial: u64, field: &str) -> Result<&'a str> {
    std::str::from_utf8(raw)
        .map_err(|_| invalid(format!("annotation record {serial} has non-UTF-8 {field}")))
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{
        self as vcf,
        header::record::value::map::info::Number,
        variant::{
            RecordBuf,
            record_buf::{
                info::field::{Value, value::Array},
                samples::sample::Value as SampleValue,
            },
        },
    };

    use super::*;
    use crate::annotate::{
        columns::{BoundColumns, ColumnSpec, SourceKind},
        matching::Matched,
        source::{AnnotationRecord, AnnotationSource, variant_record},
    };

    const BASE: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n";

    fn header(lines: &str) -> vcf::Header {
        format!("{BASE}{lines}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            .parse()
            .unwrap()
    }

    fn sample_header(lines: &str, samples: &[&str]) -> vcf::Header {
        format!(
            "{BASE}{lines}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
            samples.join("\t")
        )
        .parse()
        .unwrap()
    }

    fn record(header: &vcf::Header, line: &str) -> RecordBuf {
        let raw = vcf::Record::try_from(line.as_bytes()).unwrap();
        RecordBuf::try_from_variant_record(header, &raw).unwrap()
    }

    fn bind(target: &mut vcf::Header, source: &vcf::Header, columns: &str) -> Editor {
        let spec = ColumnSpec::parse(columns).unwrap();
        let columns = BoundColumns::bind(spec, SourceKind::Variant, target, Some(source)).unwrap();
        Editor::bind(target, Some(source), columns).unwrap()
    }

    fn annotation(header: &vcf::Header, line: &str) -> AnnotationRecord {
        variant_record(record(header, line), 1).unwrap()
    }

    fn integer_array(values: &[Option<i32>]) -> Value {
        Value::Array(Array::Integer(values.to_vec()))
    }

    fn table(
        target: &vcf::Header,
        content: &str,
        columns: &str,
    ) -> (tempfile::TempDir, AnnotationSource) {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("annotations.tsv");
        std::fs::write(&path, content).unwrap();
        let spec = ColumnSpec::parse(columns).unwrap();
        let source = AnnotationSource::open(&path, target, spec).unwrap();
        (directory, source)
    }

    #[test]
    fn copies_renamed_info_and_filter_definitions() {
        let source = header(
            "##FILTER=<ID=q10,Description=\"Low quality\">\n\
             ##INFO=<ID=OLD,Number=1,Type=Integer,Description=\"Old\">\n",
        );
        let mut target = header("");
        let _editor = bind(&mut target, &source, "FILTER,INFO/NEW:=INFO/OLD");

        assert!(target.filters().contains_key("q10"));
        assert_eq!(
            target.infos().get("NEW").unwrap().number(),
            Number::Count(1)
        );
        assert!(!target.infos().contains_key("OLD"));
    }

    #[test]
    fn selects_samples_by_name_in_source_order() {
        let source = sample_header("", &["A", "B", "C"]);
        let target = sample_header("", &["C", "B", "A"]);
        assert_eq!(
            SampleSelection::bind(&source, &target, None, None)
                .unwrap()
                .pairs(),
            &[(0, 2), (1, 1), (2, 0)]
        );
        assert_eq!(
            SampleSelection::bind(&source, &target, Some("A,C"), None)
                .unwrap()
                .pairs(),
            &[(0, 2), (2, 0)]
        );
        assert_eq!(
            SampleSelection::bind(&source, &target, Some("^B"), None)
                .unwrap()
                .pairs(),
            &[(0, 2), (2, 0)]
        );
    }

    #[test]
    fn validates_sample_selection_and_files() {
        let source = sample_header("", &["A", "B"]);
        let target = sample_header("", &["B", "A"]);
        assert!(SampleSelection::bind(&source, &target, Some("missing"), None).is_err());
        assert!(
            SampleSelection::bind(&source, &sample_header("", &["A"]), Some("B"), None).is_err()
        );

        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("samples.txt");
        std::fs::write(&path, "B\n").unwrap();
        assert_eq!(
            SampleSelection::bind(&source, &target, None, Some(&path))
                .unwrap()
                .pairs(),
            &[(1, 0)]
        );
        assert!(SampleSelection::bind(&source, &target, Some("A"), Some(&path)).is_err());
    }

    #[test]
    fn format_transfers_samples_by_name_and_remaps_gt_r_and_g() {
        let definitions = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n";
        let source_header = sample_header(definitions, &["A", "B"]);
        let mut target_header = sample_header(definitions, &["B", "A"]);
        let editor = bind(
            &mut target_header,
            &source_header,
            "FORMAT/GT,FORMAT/DP,FORMAT/AD,FORMAT/PL",
        );
        let source = annotation(
            &source_header,
            "chr1\t1\t.\tA\tG,C\t.\tPASS\t.\tGT:DP:AD:PL\t1|2:7:10,2,3:0,1,2,3,4,5\t0/1:8:11,4,5:10,11,12,13,14,15",
        );
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(2), Some(1)],
        };
        let mut target = record(
            &target_header,
            "chr1\t1\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/0\t0/0",
        );

        assert!(
            editor
                .apply_samples(&target_header, &matched, &mut target)
                .unwrap()
        );
        let a = target.samples().get(&target_header, "A").unwrap();
        assert_eq!(
            a.get("GT"),
            Some(Some(&SampleValue::Genotype("2|1".parse().unwrap())))
        );
        assert_eq!(a.get("DP"), Some(Some(&SampleValue::Integer(7))));
        assert_eq!(
            a.get("AD"),
            Some(Some(&SampleValue::from(vec![Some(10), Some(3), Some(2)])))
        );
        assert_eq!(
            a.get("PL"),
            Some(Some(&SampleValue::from(vec![
                Some(0),
                Some(3),
                Some(5),
                Some(1),
                Some(4),
                Some(2),
            ])))
        );
        let b = target.samples().get(&target_header, "B").unwrap();
        assert_eq!(b.get("DP"), Some(Some(&SampleValue::Integer(8))));
    }

    #[test]
    fn format_namespace_excludes_gt_and_initializes_new_keys() {
        let definitions = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n";
        let source_header = sample_header(definitions, &["A"]);
        let mut target_header = sample_header(definitions, &["A", "B"]);
        let editor = bind(&mut target_header, &source_header, "FORMAT");
        let source = annotation(&source_header, "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT:DP\t1/1:9");
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(1)],
        };
        let mut target = record(&target_header, "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1");

        editor
            .apply_samples(&target_header, &matched, &mut target)
            .unwrap();
        let a = target.samples().get(&target_header, "A").unwrap();
        assert_eq!(
            a.get("GT"),
            Some(Some(&SampleValue::Genotype("0/0".parse().unwrap())))
        );
        assert_eq!(a.get("DP"), Some(Some(&SampleValue::Integer(9))));
        let b = target.samples().get(&target_header, "B").unwrap();
        assert_eq!(
            b.get("GT"),
            Some(Some(&SampleValue::Genotype("0/1".parse().unwrap())))
        );
        assert_eq!(b.get("DP"), Some(None));
    }

    #[test]
    fn format_respects_selection_and_write_modes() {
        let definitions = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"X\">\n";
        let source_header = sample_header(definitions, &["A", "B"]);
        let cases = [
            ("FORMAT/X", "7", "1", Some(Some(7))),
            (".FORMAT/X", ".", "1", Some(None)),
            ("+FORMAT/X", "7", "1", Some(Some(1))),
            (".+FORMAT/X", "7", ".", Some(Some(7))),
            ("+FORMAT/X", "7", ".", Some(Some(7))),
        ];
        for (columns, source_x, target_x, expected) in cases {
            let mut target_header = sample_header(definitions, &["B", "A"]);
            let mut editor = bind(&mut target_header, &source_header, columns);
            editor.set_samples(
                SampleSelection::bind(&source_header, &target_header, Some("A"), None).unwrap(),
            );
            let source = annotation(
                &source_header,
                &format!("chr1\t1\t.\tA\tC\t.\tPASS\t.\tX\t{source_x}\t8"),
            );
            let matched = Matched {
                source: &source,
                allele_map: vec![Some(0), Some(1)],
            };
            let mut target = record(
                &target_header,
                &format!("chr1\t1\t.\tA\tC\t.\tPASS\t.\tX\t2\t{target_x}"),
            );

            editor
                .apply_samples(&target_header, &matched, &mut target)
                .unwrap();
            let a = target.samples().get(&target_header, "A").unwrap();
            let actual = a.get("X").map(|value| {
                value.map(|value| match value {
                    SampleValue::Integer(value) => *value,
                    _ => panic!("unexpected FORMAT/X value"),
                })
            });
            assert_eq!(actual, expected, "{columns}");
            let b = target.samples().get(&target_header, "B").unwrap();
            assert_eq!(
                b.get("X"),
                Some(Some(&SampleValue::Integer(2))),
                "{columns}"
            );
        }

        for columns in ["=FORMAT/X", ".=FORMAT/X"] {
            let mut target_header = sample_header(definitions, &["A"]);
            let spec = ColumnSpec::parse(columns).unwrap();
            let columns = BoundColumns::bind(
                spec,
                SourceKind::Variant,
                &target_header,
                Some(&source_header),
            )
            .unwrap();
            assert!(Editor::bind(&mut target_header, Some(&source_header), columns).is_err());
        }
    }

    #[test]
    fn format_replace_existing_does_not_create_an_absent_key() {
        let definitions = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"X\">\n";
        let source_header = sample_header(definitions, &["A"]);
        let mut target_header = sample_header(definitions, &["A"]);
        let editor = bind(&mut target_header, &source_header, "-FORMAT/X");
        let source = annotation(&source_header, "chr1\t1\t.\tA\tC\t.\tPASS\t.\tX\t7");
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(1)],
        };
        let mut target = record(&target_header, "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1");

        assert!(
            !editor
                .apply_samples(&target_header, &matched, &mut target)
                .unwrap()
        );
        assert!(!target.samples().keys().as_ref().contains("X"));
    }

    #[test]
    fn format_remaps_a_r_g_for_mixed_ploidy() {
        let definitions = "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"AO\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n";
        let source_header = sample_header(definitions, &["H", "D", "T"]);
        let mut target_header = source_header.clone();
        let editor = bind(
            &mut target_header,
            &source_header,
            "FORMAT/AO,FORMAT/AD,FORMAT/PL",
        );
        let source = annotation(
            &source_header,
            "chr1\t1\t.\tA\tG,C\t.\tPASS\t.\tAO:AD:PL\t1,2:10,1,2:0,1,2\t3,4:20,3,4:0,1,2,3,4,5\t5,6:30,5,6:0,1,2,3,4,5,6,7,8,9",
        );
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(2), Some(1)],
        };
        let mut target = record(
            &target_header,
            "chr1\t1\t.\tA\tC,G\t.\tPASS\t.\tAO:AD:PL\t.:.:.\t.:.:.\t.:.:.",
        );

        editor
            .apply_samples(&target_header, &matched, &mut target)
            .unwrap();
        let h = target.samples().get(&target_header, "H").unwrap();
        assert_eq!(
            h.get("AO"),
            Some(Some(&SampleValue::from(vec![Some(2), Some(1)])))
        );
        assert_eq!(
            h.get("AD"),
            Some(Some(&SampleValue::from(vec![Some(10), Some(2), Some(1)])))
        );
        assert_eq!(
            h.get("PL"),
            Some(Some(&SampleValue::from(vec![Some(0), Some(2), Some(1)])))
        );
        let d = target.samples().get(&target_header, "D").unwrap();
        assert_eq!(
            d.get("PL"),
            Some(Some(&SampleValue::from(vec![
                Some(0),
                Some(3),
                Some(5),
                Some(1),
                Some(4),
                Some(2)
            ])))
        );
        let t = target.samples().get(&target_header, "T").unwrap();
        assert_eq!(
            t.get("PL"),
            Some(Some(&SampleValue::from(vec![
                Some(0),
                Some(4),
                Some(7),
                Some(9),
                Some(1),
                Some(5),
                Some(8),
                Some(2),
                Some(6),
                Some(3)
            ])))
        );
    }

    #[test]
    fn format_preserves_gt_separators_and_rejects_unmapped_alleles() {
        let definitions = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n";
        let source_header = sample_header(definitions, &["A"]);
        let mut target_header = source_header.clone();
        let editor = bind(&mut target_header, &source_header, "FORMAT/GT");
        let source = annotation(&source_header, "chr1\t1\t.\tA\tG,C\t.\tPASS\t.\tGT\t1|2/1");
        let mut target = record(&target_header, "chr1\t1\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/0");
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(2), Some(1)],
        };
        editor
            .apply_samples(&target_header, &matched, &mut target)
            .unwrap();
        assert_eq!(
            target.samples().get(&target_header, "A").unwrap().get("GT"),
            Some(Some(&SampleValue::Genotype("2|1/2".parse().unwrap())))
        );

        let impossible = Matched {
            source: &source,
            allele_map: vec![Some(0), None, Some(1)],
        };
        assert!(
            editor
                .apply_samples(&target_header, &impossible, &mut target)
                .is_err()
        );
    }

    #[test]
    fn format_rejects_invalid_schema_and_cardinality() {
        let source_header = sample_header(
            "##FORMAT=<ID=X,Number=2,Type=Integer,Description=\"X\">\n",
            &["A"],
        );
        let mut target_header = source_header.clone();
        let editor = bind(&mut target_header, &source_header, "FORMAT/X");
        let source = annotation(&source_header, "chr1\t1\t.\tA\tC\t.\tPASS\t.\tX\t7");
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(1)],
        };
        let mut target = record(&target_header, "chr1\t1\t.\tA\tC\t.\tPASS\t.\tX\t1,2");
        assert!(
            editor
                .apply_samples(&target_header, &matched, &mut target)
                .is_err()
        );

        let source_header = sample_header(
            "##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"X\">\n",
            &["A"],
        );
        let mut target_header = sample_header(
            "##FORMAT=<ID=X,Number=1,Type=String,Description=\"X\">\n",
            &["A"],
        );
        let spec = ColumnSpec::parse("FORMAT/X").unwrap();
        let columns = BoundColumns::bind(
            spec,
            SourceKind::Variant,
            &target_header,
            Some(&source_header),
        )
        .unwrap();
        assert!(Editor::bind(&mut target_header, Some(&source_header), columns).is_err());
    }

    #[test]
    fn rejects_incompatible_existing_info_schema() {
        let source = header("##INFO=<ID=X,Number=1,Type=Integer,Description=\"X\">\n");
        let mut target = header("##INFO=<ID=X,Number=1,Type=String,Description=\"X\">\n");
        let spec = ColumnSpec::parse("INFO/X").unwrap();
        let columns =
            BoundColumns::bind(spec, SourceKind::Variant, &target, Some(&source)).unwrap();
        assert!(Editor::bind(&mut target, Some(&source), columns).is_err());

        let mut target = header("");
        let spec = ColumnSpec::parse("ID:=INFO/X").unwrap();
        let columns =
            BoundColumns::bind(spec, SourceKind::Variant, &target, Some(&source)).unwrap();
        assert!(Editor::bind(&mut target, Some(&source), columns).is_err());
    }

    #[test]
    fn header_preparation_is_atomic_and_rejects_effective_duplicates() {
        let source = header(
            "##INFO=<ID=A,Number=1,Type=Integer,Description=\"A\">\n\
             ##INFO=<ID=B,Number=1,Type=Integer,Description=\"B\">\n",
        );
        let mut target = header("##INFO=<ID=B,Number=1,Type=String,Description=\"B\">\n");
        let before = target.clone();
        let spec = ColumnSpec::parse("INFO").unwrap();
        let columns =
            BoundColumns::bind(spec, SourceKind::Variant, &target, Some(&source)).unwrap();
        assert!(Editor::bind(&mut target, Some(&source), columns).is_err());
        assert_eq!(target, before);

        let mut target = header("");
        let spec = ColumnSpec::parse("INFO,INFO/A").unwrap();
        let columns =
            BoundColumns::bind(spec, SourceKind::Variant, &target, Some(&source)).unwrap();
        assert!(Editor::bind(&mut target, Some(&source), columns).is_err());
        assert!(!target.infos().contains_key("A"));
    }

    #[test]
    fn transfers_fixed_and_typed_info_fields() {
        let definitions = "##FILTER=<ID=q10,Description=\"Low quality\">\n\
##INFO=<ID=FLAG,Number=0,Type=Flag,Description=\"Flag\">\n\
##INFO=<ID=I,Number=1,Type=Integer,Description=\"Integer\">\n\
##INFO=<ID=F,Number=1,Type=Float,Description=\"Float\">\n\
##INFO=<ID=C,Number=1,Type=Character,Description=\"Character\">\n\
##INFO=<ID=S,Number=1,Type=String,Description=\"String\">\n";
        let source_header = header(definitions);
        let mut target_header = header(definitions);
        let editor = bind(
            &mut target_header,
            &source_header,
            "ID,QUAL,FILTER,INFO/FLAG,INFO/I,INFO/F,INFO/C,INFO/S",
        );
        let source = annotation(
            &source_header,
            "chr1\t1\ts1;s2\tA\tC\t42\tq10\tFLAG;I=7;F=1.5;C=Z;S=text",
        );
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(1)],
        };
        let mut target = record(&target_header, "chr1\t1\told\tA\tC\t5\tPASS\tI=1;S=old");

        assert!(
            editor
                .apply_info(&target_header, &matched, &mut target)
                .unwrap()
        );
        assert_eq!(
            target.ids().as_ref().iter().cloned().collect::<Vec<_>>(),
            ["s1", "s2"]
        );
        assert_eq!(target.quality_score(), Some(42.0));
        assert_eq!(
            target
                .filters()
                .as_ref()
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            ["q10"]
        );
        assert_eq!(target.info().get("FLAG"), Some(Some(&Value::Flag)));
        assert_eq!(target.info().get("I"), Some(Some(&Value::Integer(7))));
        assert_eq!(target.info().get("F"), Some(Some(&Value::Float(1.5))));
        assert_eq!(target.info().get("C"), Some(Some(&Value::Character('Z'))));
        assert_eq!(
            target.info().get("S"),
            Some(Some(&Value::String("text".into())))
        );
    }

    #[test]
    fn parses_tabular_fixed_fields_and_deduplicates_sets() {
        let mut target_header = header("##FILTER=<ID=q10,Description=\"Low quality\">\n");
        let (_directory, mut source) = table(
            &target_header,
            "chr1\t1\ts1;s1;s2\t42\tq10;q10\n",
            "CHROM,POS,ID,QUAL,FILTER",
        );
        let editor = Editor::bind(&mut target_header, None, source.columns().clone()).unwrap();
        let annotation = source.read_record().unwrap().unwrap();
        let matched = Matched {
            source: &annotation,
            allele_map: Vec::new(),
        };
        let mut target = record(&target_header, "chr1\t1\told\tA\tC\t5\tPASS\t.");
        editor
            .apply_info(&target_header, &matched, &mut target)
            .unwrap();

        assert_eq!(
            target.ids().as_ref().iter().cloned().collect::<Vec<_>>(),
            ["s1", "s2"]
        );
        assert_eq!(target.quality_score(), Some(42.0));
        assert_eq!(
            target
                .filters()
                .as_ref()
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            ["q10"]
        );
    }

    #[test]
    fn applies_replace_add_missing_and_existing_modes() {
        let source_header = header("##INFO=<ID=X,Number=1,Type=Integer,Description=\"X\">\n");
        let cases = [
            ("INFO/X", "X=7", "X=1", Some(Some(&Value::Integer(7)))),
            (".INFO/X", "X=.", "X=1", Some(None)),
            ("+INFO/X", "X=7", "X=1", Some(Some(&Value::Integer(1)))),
            (".+INFO/X", "X=7", "X=.", Some(Some(&Value::Integer(7)))),
            ("-INFO/X", "X=7", ".", None),
            ("-INFO/X", "X=7", "X=1", Some(Some(&Value::Integer(7)))),
        ];
        for (columns, source_info, target_info, expected) in cases {
            let mut target_header = source_header.clone();
            let editor = bind(&mut target_header, &source_header, columns);
            let source = annotation(
                &source_header,
                &format!("chr1\t1\t.\tA\tC\t.\tPASS\t{source_info}"),
            );
            let matched = Matched {
                source: &source,
                allele_map: vec![Some(0), Some(1)],
            };
            let mut target = record(
                &target_header,
                &format!("chr1\t1\t.\tA\tC\t.\tPASS\t{target_info}"),
            );
            editor
                .apply_info(&target_header, &matched, &mut target)
                .unwrap();
            assert_eq!(target.info().get("X"), expected, "{columns}");
        }
    }

    #[test]
    fn transfers_all_info_fields_and_missing_values() {
        let source_header = header(
            "##INFO=<ID=A,Number=1,Type=Integer,Description=\"A\">\n\
             ##INFO=<ID=B,Number=1,Type=String,Description=\"B\">\n",
        );
        let mut target_header = header("");
        let editor = bind(&mut target_header, &source_header, "INFO");
        let source = annotation(&source_header, "chr1\t1\t.\tA\tC\t.\tPASS\tA=7;B=.");
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(1)],
        };
        let mut target = record(&target_header, "chr1\t1\t.\tA\tC\t.\tPASS\t.");

        editor
            .apply_info(&target_header, &matched, &mut target)
            .unwrap();
        assert_eq!(target.info().get("A"), Some(Some(&Value::Integer(7))));
        assert_eq!(target.info().get("B"), Some(None));
        assert!(target_header.infos().contains_key("A"));
        assert!(target_header.infos().contains_key("B"));
    }

    #[test]
    fn parses_tabular_values_and_applies_append_modes() {
        let mut target_header = header(
            "##INFO=<ID=X,Number=1,Type=Integer,Description=\"X\">\n\
             ##INFO=<ID=FLAG,Number=0,Type=Flag,Description=\"Flag\">\n",
        );
        let (_directory, mut source) = table(
            &target_header,
            "chr1\t1\t7,8\t1\n",
            "CHROM,POS,=INFO/X,INFO/FLAG",
        );
        let editor = Editor::bind(&mut target_header, None, source.columns().clone()).unwrap();
        assert_eq!(
            target_header.infos().get("X").unwrap().number(),
            Number::Unknown
        );
        let annotation = source.read_record().unwrap().unwrap();
        let matched = Matched {
            source: &annotation,
            allele_map: Vec::new(),
        };
        let mut target = record(&target_header, "chr1\t1\t.\tA\tC\t.\tPASS\tX=1,7");
        editor
            .apply_info(&target_header, &matched, &mut target)
            .unwrap();
        assert_eq!(
            target.info().get("X"),
            Some(Some(&integer_array(&[Some(1), Some(7), Some(7), Some(8)])))
        );

        for (raw, expected) in [
            ("7,8", integer_array(&[None, Some(7), Some(8)])),
            (".", integer_array(&[None, None])),
        ] {
            let (_directory, mut source) = table(
                &target_header,
                &format!("chr1\t1\t{raw}\n"),
                "CHROM,POS,.=INFO/X",
            );
            let editor = Editor::bind(&mut target_header, None, source.columns().clone()).unwrap();
            let annotation = source.read_record().unwrap().unwrap();
            let matched = Matched {
                source: &annotation,
                allele_map: Vec::new(),
            };
            let mut missing = record(&target_header, "chr1\t1\t.\tA\tC\t.\tPASS\tX=.");
            editor
                .apply_info(&target_header, &matched, &mut missing)
                .unwrap();
            assert_eq!(missing.info().get("X"), Some(Some(&expected)), "{raw}");
        }
        assert_eq!(target.info().get("FLAG"), Some(Some(&Value::Flag)));

        let (_directory, mut source) = table(&target_header, "chr1\t1\t.\n", "CHROM,POS,.=INFO/X");
        let editor = Editor::bind(&mut target_header, None, source.columns().clone()).unwrap();
        let annotation = source.read_record().unwrap().unwrap();
        let matched = Matched {
            source: &annotation,
            allele_map: Vec::new(),
        };
        editor
            .apply_info(&target_header, &matched, &mut target)
            .unwrap();
        assert_eq!(
            target.info().get("X"),
            Some(Some(&integer_array(&[
                Some(1),
                Some(7),
                Some(7),
                Some(8),
                None,
            ])))
        );
    }

    #[test]
    fn tabular_missing_values_require_dot_write_modes() {
        let base = header("##INFO=<ID=X,Number=1,Type=Integer,Description=\"X\">\n");
        for (columns, expected) in [
            ("CHROM,POS,INFO/X", Some(Some(&Value::Integer(1)))),
            ("CHROM,POS,.INFO/X", Some(None)),
        ] {
            let mut target_header = base.clone();
            let (_directory, mut source) = table(&target_header, "chr1\t1\t.\n", columns);
            let editor = Editor::bind(&mut target_header, None, source.columns().clone()).unwrap();
            let annotation = source.read_record().unwrap().unwrap();
            let matched = Matched {
                source: &annotation,
                allele_map: Vec::new(),
            };
            let mut target = record(&target_header, "chr1\t1\t.\tA\tC\t.\tPASS\tX=1");
            editor
                .apply_info(&target_header, &matched, &mut target)
                .unwrap();
            assert_eq!(target.info().get("X"), expected, "{columns}");
        }
    }

    #[test]
    fn tabular_flags_follow_bcftools_three_state_semantics() {
        let base = header("##INFO=<ID=F,Number=0,Type=Flag,Description=\"F\">\n");
        let cases = [
            ("0", "CHROM,POS,INFO/F", "F", None),
            (".", "CHROM,POS,.INFO/F", "F", Some(Some(&Value::Flag))),
            ("1", "CHROM,POS,+INFO/F", ".", Some(Some(&Value::Flag))),
        ];
        for (raw, columns, target_info, expected) in cases {
            let mut target_header = base.clone();
            let (_directory, mut source) =
                table(&target_header, &format!("chr1\t1\t{raw}\n"), columns);
            let editor = Editor::bind(&mut target_header, None, source.columns().clone()).unwrap();
            let annotation = source.read_record().unwrap().unwrap();
            let matched = Matched {
                source: &annotation,
                allele_map: Vec::new(),
            };
            let mut target = record(
                &target_header,
                &format!("chr1\t1\t.\tA\tC\t.\tPASS\t{target_info}"),
            );
            editor
                .apply_info(&target_header, &matched, &mut target)
                .unwrap();
            assert_eq!(target.info().get("F"), expected, "{raw}");
        }

        let mut target_header = base;
        let (_directory, source) = table(&target_header, "chr1\t1\t1\n", "CHROM,POS,=INFO/F");
        assert!(Editor::bind(&mut target_header, None, source.columns().clone()).is_err());
    }

    #[test]
    fn remaps_number_a_r_and_g_for_reordered_alleles() {
        let source_header = header("");
        let source = annotation(&source_header, "chr1\t1\t.\tA\tG,C\t.\tPASS\t.");
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(2), Some(1)],
        };

        assert_eq!(
            remap_info(
                Number::A,
                integer_array(&[Some(2), Some(3)]),
                &matched,
                3,
                "INFO/A"
            )
            .unwrap(),
            integer_array(&[Some(3), Some(2)])
        );
        assert_eq!(
            remap_info(
                Number::R,
                integer_array(&[Some(10), Some(2), Some(3)]),
                &matched,
                3,
                "INFO/R",
            )
            .unwrap(),
            integer_array(&[Some(10), Some(3), Some(2)])
        );
        assert_eq!(
            remap_info(
                Number::G,
                integer_array(&[Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)]),
                &matched,
                3,
                "INFO/G",
            )
            .unwrap(),
            integer_array(&[Some(0), Some(3), Some(5), Some(1), Some(4), Some(2)])
        );
    }

    #[test]
    fn editor_remaps_typed_allele_fields() {
        let definitions = "##INFO=<ID=A,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=R,Number=R,Type=Integer,Description=\"R\">\n\
##INFO=<ID=G,Number=G,Type=Integer,Description=\"G\">\n";
        let source_header = header(definitions);
        let mut target_header = header(definitions);
        let editor = bind(&mut target_header, &source_header, "INFO/A,INFO/R,INFO/G");
        let source = annotation(
            &source_header,
            "chr1\t1\t.\tA\tG,C\t.\tPASS\tA=2,3;R=10,2,3;G=0,1,2,3,4,5",
        );
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(2), Some(1)],
        };
        let mut target = record(&target_header, "chr1\t1\t.\tA\tC,G\t.\tPASS\t.");
        editor
            .apply_info(&target_header, &matched, &mut target)
            .unwrap();

        assert_eq!(
            target.info().get("A"),
            Some(Some(&integer_array(&[Some(3), Some(2)])))
        );
        assert_eq!(
            target.info().get("R"),
            Some(Some(&integer_array(&[Some(10), Some(3), Some(2)])))
        );
        assert_eq!(
            target.info().get("G"),
            Some(Some(&integer_array(&[
                Some(0),
                Some(3),
                Some(5),
                Some(1),
                Some(4),
                Some(2),
            ])))
        );
    }

    #[test]
    fn remaps_subset_extended_and_multiple_ploidies() {
        let source_header = header("");
        let source = annotation(&source_header, "chr1\t1\t.\tA\tG,C\t.\tPASS\t.");
        let subset = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(1), None],
        };
        assert_eq!(
            remap_info(
                Number::A,
                integer_array(&[Some(2), Some(3)]),
                &subset,
                2,
                "INFO/A"
            )
            .unwrap(),
            integer_array(&[Some(2)])
        );
        assert_eq!(
            remap_info(
                Number::G,
                integer_array(&[Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)]),
                &subset,
                2,
                "INFO/G",
            )
            .unwrap(),
            integer_array(&[Some(0), Some(1), Some(2)])
        );

        let source = annotation(&source_header, "chr1\t1\t.\tA\tG\t.\tPASS\t.");
        let extended = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(2)],
        };
        assert_eq!(
            remap_info(
                Number::R,
                integer_array(&[Some(10), Some(2)]),
                &extended,
                3,
                "INFO/R",
            )
            .unwrap(),
            integer_array(&[Some(10), None, Some(2)])
        );

        let identity = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(1)],
        };
        for values in [2, 3, 4] {
            let input = integer_array(&(0..values).map(Some).collect::<Vec<_>>());
            assert_eq!(
                remap_info(Number::G, input.clone(), &identity, 2, "INFO/G").unwrap(),
                input
            );
        }
    }

    #[test]
    fn remaps_float_character_and_string_arrays_without_coercion() {
        let source_header = header("");
        let source = annotation(&source_header, "chr1\t1\t.\tA\tG,C\t.\tPASS\t.");
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(2), Some(1)],
        };
        let values = [
            (
                Value::Array(Array::Float(vec![Some(1.0), Some(2.0), Some(3.0)])),
                Value::Array(Array::Float(vec![Some(1.0), Some(3.0), Some(2.0)])),
            ),
            (
                Value::Array(Array::Character(vec![Some('R'), Some('G'), Some('C')])),
                Value::Array(Array::Character(vec![Some('R'), Some('C'), Some('G')])),
            ),
            (
                Value::Array(Array::String(vec![
                    Some("R".into()),
                    Some("G".into()),
                    Some("C".into()),
                ])),
                Value::Array(Array::String(vec![
                    Some("R".into()),
                    Some("C".into()),
                    Some("G".into()),
                ])),
            ),
        ];
        for (source, expected) in values {
            assert_eq!(
                remap_info(Number::R, source, &matched, 3, "INFO/X").unwrap(),
                expected
            );
        }
    }

    #[test]
    fn rejects_invalid_allele_cardinality_and_value_types() {
        let source_header = header("");
        let source = annotation(&source_header, "chr1\t1\t.\tA\tG,C\t.\tPASS\t.");
        let matched = Matched {
            source: &source,
            allele_map: vec![Some(0), Some(2), Some(1)],
        };
        assert!(remap_info(Number::A, Value::Integer(1), &matched, 3, "INFO/A").is_err());
        assert!(remap_info(Number::R, integer_array(&[Some(1)]), &matched, 3, "INFO/R").is_err());
        assert!(
            remap_info(
                Number::G,
                integer_array(&[Some(1), Some(2), Some(3), Some(4)]),
                &matched,
                3,
                "INFO/G",
            )
            .is_err()
        );
    }
}
