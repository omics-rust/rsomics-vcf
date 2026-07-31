mod reader;
mod record;
mod text;
mod value;

pub(crate) use reader::Reader;
pub(crate) use record::reformat_record;
pub(crate) use text::HeaderTypes;
pub(crate) use value::{ValueType, write_f32};

pub(crate) fn trim_line_ending(line: &mut Vec<u8>) {
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
}
