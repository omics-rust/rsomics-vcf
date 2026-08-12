use noodles_vcf::variant::{
    RecordBuf,
    record_buf::info::field::{Value as InfoValue, value::Array as InfoArray},
};

use crate::variant_type;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Policy {
    Snps,
    Indels,
    Both,
    All,
    Exact,
}

pub(super) struct State {
    policy: Option<Policy>,
    coordinate: Option<(usize, usize)>,
    types: u32,
    alleles: Vec<Alleles>,
}

struct Alleles {
    reference: String,
    alternates: Vec<String>,
}

impl State {
    pub(super) fn new(policy: Option<Policy>) -> Self {
        Self {
            policy,
            coordinate: None,
            types: 0,
            alleles: Vec::new(),
        }
    }

    pub(super) fn remove(&mut self, coordinate: (usize, usize), record: &RecordBuf) -> bool {
        let Some(policy) = self.policy else {
            return false;
        };
        if self.coordinate != Some(coordinate) {
            self.coordinate = Some(coordinate);
            self.types = 0;
            self.alleles.clear();
        }

        let types = variant_type::record_mask(record);
        let duplicate = match policy {
            Policy::All => self.types != 0,
            Policy::Snps => same_type(types, self.types, variant_type::SNP | variant_type::MNP),
            Policy::Indels => same_type(types, self.types, variant_type::INDEL),
            Policy::Both => {
                same_type(types, self.types, variant_type::SNP | variant_type::MNP)
                    || same_type(types, self.types, variant_type::INDEL)
            }
            Policy::Exact => {
                let alleles = Alleles::from_record(record);
                if self
                    .alleles
                    .iter()
                    .any(|previous| previous.matches(&alleles))
                {
                    true
                } else {
                    self.alleles.push(alleles);
                    false
                }
            }
        };
        if !duplicate {
            self.types |= types;
        }
        duplicate
    }
}

impl Alleles {
    fn from_record(record: &RecordBuf) -> Self {
        let raw = record.alternate_bases().as_ref();
        let alternates = raw
            .iter()
            .enumerate()
            .map(|(index, alternate)| {
                let mut alternate = alternate.to_ascii_uppercase();
                if raw.len() == 1
                    && alternate.starts_with('<')
                    && let Some(length) = svlen(record, index)
                {
                    alternate.push('.');
                    alternate.push_str(&length.to_string());
                }
                alternate
            })
            .collect();
        Self {
            reference: record.reference_bases().to_ascii_uppercase(),
            alternates,
        }
    }

    fn matches(&self, other: &Self) -> bool {
        self.reference == other.reference
            && self.alternates.len() == other.alternates.len()
            && if self.alternates.len() == 1 {
                self.alternates == other.alternates
            } else {
                other
                    .alternates
                    .iter()
                    .all(|alternate| self.alternates.contains(alternate))
            }
    }
}

fn same_type(current: u32, previous: u32, mask: u32) -> bool {
    current & mask != 0 && previous & mask != 0
}

fn svlen(record: &RecordBuf, index: usize) -> Option<i32> {
    match record.info().get("SVLEN")? {
        Some(InfoValue::Integer(value)) if index == 0 => Some(*value),
        Some(InfoValue::Array(InfoArray::Integer(values))) => values.get(index).copied().flatten(),
        _ => None,
    }
}
