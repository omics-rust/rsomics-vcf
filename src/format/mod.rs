mod reader;
mod record;
mod text;
mod value;
mod writer;

use serde::Serialize;

pub(crate) use reader::{Reader, RecordScratch};
pub(crate) use record::reformat_record;
pub(crate) use text::HeaderTypes;
pub(crate) use value::{ValueType, write_f32};
pub(crate) use writer::{ParallelWriter, VariantWriter, Writer};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum OutputFormat {
    #[default]
    Vcf,
    VcfBgzf,
    Bcf,
    BcfRaw,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum HeaderMode {
    #[default]
    Full,
    HeaderOnly,
    None,
}

pub(crate) fn trim_line_ending(line: &mut Vec<u8>) {
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
}
