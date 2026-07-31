use noodles_vcf::variant::RecordBuf;
use rsomics_common::{Result, RsomicsError};

use super::{IdSelection, Options};
use crate::variant_type;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TypeSelection {
    mask: u32,
    include: bool,
}

impl TypeSelection {
    pub fn include(values: &str) -> Result<Self> {
        Self::parse(values, true)
    }

    pub fn exclude(values: &str) -> Result<Self> {
        Self::parse(values, false)
    }

    fn parse(values: &str, include: bool) -> Result<Self> {
        let mask = variant_type::parse_mask(values).map_err(RsomicsError::InvalidInput)?;
        Ok(Self { mask, include })
    }
}

pub(super) fn keep(record: &RecordBuf, options: &Options) -> bool {
    let alleles = record.alternate_bases().as_ref();
    let allele_count = alleles.len() + 1;
    if options
        .min_alleles
        .is_some_and(|minimum| allele_count < minimum)
        || options
            .max_alleles
            .is_some_and(|maximum| allele_count > maximum)
    {
        return false;
    }

    if let Some(selection) = options.ids {
        let record_known = !record.ids().as_ref().is_empty();
        if record_known != (selection == IdSelection::Known) {
            return false;
        }
    }

    if !options.filters.is_empty() {
        let values = record.filters().as_ref();
        let matched = options.filters.iter().any(|filter| {
            (filter == "PASS" && record.filters().is_pass())
                || (filter == "." && values.is_empty())
                || values.contains(filter)
        });
        if !matched {
            return false;
        }
    }

    if let Some(types) = options.types {
        let joined = alleles.join(",");
        let matched = variant_type::mask(record.reference_bases().as_bytes(), joined.as_bytes())
            & types.mask
            != 0;
        if matched != types.include {
            return false;
        }
    }

    true
}
