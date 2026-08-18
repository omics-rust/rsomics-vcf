use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use rsomics_common::{Context, Result, RsomicsError};

pub(super) struct Fai {
    entries: Vec<FaiEntry>,
}

struct FaiEntry {
    name: String,
    length: u64,
}

impl Fai {
    pub(super) fn read(path: &Path) -> Result<Self> {
        let source =
            fs::read(path).rs_with_context(|| format!("reading FASTA index {}", path.display()))?;
        Self::parse(&source)
    }

    fn parse(source: &[u8]) -> Result<Self> {
        let source = std::str::from_utf8(source)
            .map_err(|error| invalid(format!("FASTA index is not valid UTF-8: {error}")))?;
        let mut entries = Vec::new();
        let mut names = HashSet::new();
        for (index, line) in source.lines().enumerate() {
            if line.is_empty() {
                return Err(invalid(format!("FASTA index line {} is empty", index + 1)));
            }
            let mut fields = line.split('\t');
            let name = fields.next().unwrap_or_default();
            let length = fields.next().ok_or_else(|| {
                invalid(format!(
                    "FASTA index line {} has fewer than two columns",
                    index + 1
                ))
            })?;
            if name.is_empty() {
                return Err(invalid(format!(
                    "FASTA index line {} has an empty contig name",
                    index + 1
                )));
            }
            if name.chars().any(char::is_control) {
                return Err(invalid(format!(
                    "FASTA index line {} has a control character in the contig name",
                    index + 1
                )));
            }
            if !names.insert(name) {
                return Err(invalid(format!("duplicate FASTA index contig: {name}")));
            }
            let length = length.parse().map_err(|error| {
                invalid(format!(
                    "FASTA index line {} has invalid length {length:?}: {error}",
                    index + 1
                ))
            })?;
            entries.push(FaiEntry {
                name: name.to_owned(),
                length,
            });
        }
        if entries.is_empty() {
            return Err(invalid("FASTA index contains no contigs"));
        }
        Ok(Self { entries })
    }
}

pub(super) fn rewrite_contigs(metadata: &[String], fai: &Fai) -> Result<Vec<String>> {
    let lengths: HashMap<_, _> = fai
        .entries
        .iter()
        .map(|entry| (entry.name.as_str(), entry.length))
        .collect();
    let mut output = Vec::with_capacity(metadata.len() + fai.entries.len());
    let mut header_names = HashSet::new();
    let mut retained = HashSet::new();
    for line in metadata {
        let Some(contig) = ContigLine::parse(line)? else {
            output.push(line.clone());
            continue;
        };
        if !header_names.insert(contig.id) {
            return Err(invalid(format!(
                "duplicate VCF header contig: {}",
                contig.id
            )));
        }
        if let Some(length) = lengths.get(contig.id) {
            output.push(contig.with_length(*length));
            retained.insert(contig.id);
        }
    }
    for entry in &fai.entries {
        if !retained.contains(entry.name.as_str()) {
            output.push(format!(
                "##contig=<ID={},length={}>",
                entry.name, entry.length
            ));
        }
    }
    Ok(output)
}

struct ContigLine<'a> {
    id: &'a str,
    fields: Vec<&'a str>,
}

impl<'a> ContigLine<'a> {
    fn parse(line: &'a str) -> Result<Option<Self>> {
        let Some(body) = line.strip_prefix("##contig=<") else {
            return Ok(None);
        };
        let body = body
            .strip_suffix('>')
            .ok_or_else(|| invalid("malformed VCF contig metadata"))?;
        let fields = structured_fields(body)?;
        let mut id = None;
        let mut length = false;
        for field in &fields {
            if let Some(value) = field.strip_prefix("ID=") {
                if id.replace(value).is_some() {
                    return Err(invalid("VCF contig metadata contains more than one ID"));
                }
            } else if field.starts_with("length=") && std::mem::replace(&mut length, true) {
                return Err(invalid("VCF contig metadata contains more than one length"));
            }
        }
        let id = id.ok_or_else(|| invalid("VCF contig metadata is missing ID"))?;
        Ok(Some(Self { id, fields }))
    }

    fn with_length(&self, length: u64) -> String {
        let mut fields = Vec::with_capacity(self.fields.len() + 1);
        let mut replaced = false;
        for field in &self.fields {
            if field.starts_with("length=") {
                fields.push(format!("length={length}"));
                replaced = true;
            } else {
                fields.push((*field).to_owned());
            }
        }
        if !replaced {
            fields.push(format!("length={length}"));
        }
        format!("##contig=<{}>", fields.join(","))
    }
}

fn structured_fields(body: &str) -> Result<Vec<&str>> {
    let mut fields = Vec::new();
    let mut start = 0;
    let mut quoted = false;
    let mut escaped = false;
    for (index, character) in body.char_indices() {
        if escaped {
            escaped = false;
        } else if quoted && character == '\\' {
            escaped = true;
        } else if character == '"' {
            quoted = !quoted;
        } else if character == ',' && !quoted {
            fields.push(&body[start..index]);
            start = index + 1;
        }
    }
    if quoted || escaped {
        return Err(invalid("unterminated quoted VCF contig metadata"));
    }
    fields.push(&body[start..]);
    if fields.iter().any(|field| field.is_empty()) {
        return Err(invalid("VCF contig metadata contains an empty field"));
    }
    Ok(fields)
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use crate::reheader::header::HeaderText;

    use super::Fai;

    fn fixture(source: &[u8]) -> tempfile::NamedTempFile {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(source).unwrap();
        file
    }

    #[test]
    fn parses_first_two_columns_and_u64_boundaries() {
        let fai =
            Fai::parse(format!("zero\t0\t1\t2\t3\nmax\t{}\textra\n", u64::MAX).as_bytes()).unwrap();
        assert_eq!(fai.entries[0].name, "zero");
        assert_eq!(fai.entries[0].length, 0);
        assert_eq!(fai.entries[1].name, "max");
        assert_eq!(fai.entries[1].length, u64::MAX);
    }

    #[test]
    fn rejects_empty_duplicate_and_malformed_rows() {
        for source in [
            b"".as_slice(),
            b"\n",
            b"chr1\n",
            b"\t1\n",
            b"chr1\t\n",
            b"chr1\tx\n",
            b"chr1\t1\nchr1\t2\n",
            b"chr1\t1\n\n",
            &[0xff, b'\t', b'1', b'\n'],
        ] {
            assert!(Fai::parse(source).is_err(), "{source:?}");
        }
    }

    #[test]
    fn replaces_contigs_without_losing_other_metadata() {
        let mut header = HeaderText::parse(
            b"##fileformat=VCFv4.3\n\
##source=before\n\
##contig=<ID=chr1,length=100,assembly=old>\n\
##contig=<ID=chr2,assembly=test,Description=\"a,b\">\n\
##source=after\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        )
        .unwrap();
        let fai = fixture(b"chr2\t250\t0\t0\t0\nchr3\t300\t0\t0\t0\n");
        header.apply_fai(fai.path()).unwrap();
        assert_eq!(
            header.render(),
            b"##fileformat=VCFv4.3\n\
##source=before\n\
##contig=<ID=chr2,assembly=test,Description=\"a,b\",length=250>\n\
##source=after\n\
##contig=<ID=chr3,length=300>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        );
    }

    #[test]
    fn updates_length_in_place_and_revalidates_the_header() {
        let mut header = HeaderText::parse(
            b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100,assembly=test>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        )
        .unwrap();
        let fai = fixture(b"chr1\t18446744073709551615\n");
        header.apply_fai(fai.path()).unwrap();
        assert_eq!(
            header.render(),
            b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=18446744073709551615,assembly=test>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        );
    }

    #[test]
    fn duplicate_header_contigs_fail_before_synchronization() {
        assert!(
            HeaderText::parse(
                b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1>\n\
##contig=<ID=chr1,length=2>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            )
            .is_err()
        );
    }
}
