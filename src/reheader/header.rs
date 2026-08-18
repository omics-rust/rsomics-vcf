use std::collections::HashSet;
use std::fs;
use std::path::Path;

use noodles_vcf as vcf;
use rsomics_common::{Context, Result, RsomicsError};

use super::fai::{Fai, rewrite_contigs};
use super::samples::{SampleEdit, SampleSource};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct HeaderText {
    metadata: Vec<String>,
    columns: Vec<String>,
}

impl HeaderText {
    pub(super) fn parse(raw: &[u8]) -> Result<Self> {
        let text = std::str::from_utf8(raw)
            .map_err(|error| invalid(format!("VCF header is not valid UTF-8: {error}")))?;
        let normalized = text.replace("\r\n", "\n");
        if normalized.contains('\r') {
            return Err(invalid(
                "VCF header contains a carriage return without a line feed",
            ));
        }

        let content = normalized.strip_suffix('\n').unwrap_or(&normalized);
        let mut lines = content.split('\n');
        let fileformat = lines.next().ok_or_else(|| invalid("VCF header is empty"))?;
        if !fileformat.starts_with("##fileformat=") {
            return Err(invalid("VCF header must start with ##fileformat"));
        }
        validate_line(fileformat)?;

        let mut metadata = vec![fileformat.to_owned()];
        let mut columns = None;
        for line in lines {
            if line.is_empty() {
                return Err(invalid("VCF header contains a blank line"));
            }
            validate_line(line)?;
            if columns.is_some() {
                return Err(invalid("VCF header contains data after the #CHROM line"));
            }
            if line.starts_with("##fileformat=") {
                return Err(invalid(
                    "VCF header contains more than one ##fileformat line",
                ));
            }
            if line.starts_with("##") {
                metadata.push(line.to_owned());
            } else if line.starts_with("#CHROM") {
                columns = Some(parse_columns(line)?);
            } else {
                return Err(invalid(
                    "VCF header contains a line that is not metadata or #CHROM",
                ));
            }
        }

        let header = Self {
            metadata,
            columns: columns.ok_or_else(|| invalid("VCF header is missing the #CHROM line"))?,
        };
        header.parse_noodles()?;
        Ok(header)
    }

    pub(super) fn replace_from(path: &Path, expected_samples: usize) -> Result<Self> {
        let raw = fs::read(path)
            .rs_with_context(|| format!("reading replacement header {}", path.display()))?;
        let header = Self::parse(&raw)?;
        let actual = header.sample_names().len();
        if actual != expected_samples {
            return Err(invalid(format!(
                "replacement header must contain {expected_samples} samples, found {actual}"
            )));
        }
        Ok(header)
    }

    pub(super) fn render(&self) -> Vec<u8> {
        let mut output = self.metadata.join("\n");
        output.push('\n');
        output.push_str(&self.columns.join("\t"));
        output.push('\n');
        output.into_bytes()
    }

    pub(super) fn sample_names(&self) -> &[String] {
        if self.columns.len() > 9 {
            &self.columns[9..]
        } else {
            &[]
        }
    }

    pub(super) fn contig_count(&self) -> usize {
        self.metadata
            .iter()
            .filter(|line| line.starts_with("##contig=<"))
            .count()
    }

    pub(super) fn parse_noodles(&self) -> Result<vcf::Header> {
        let text = String::from_utf8(self.render()).expect("rendered header is UTF-8");
        text.parse()
            .map_err(|error| invalid(format!("parsing VCF header: {error}")))
    }

    pub(super) fn apply_samples(&mut self, source: &SampleSource) -> Result<()> {
        let edit = SampleEdit::read(source)?;
        let names = edit.apply(self.sample_names())?;
        self.set_samples(names)
    }

    fn set_samples(&mut self, names: Vec<String>) -> Result<()> {
        let mut columns = self.columns[..9].to_vec();
        columns.extend(names);
        let edited = Self {
            metadata: self.metadata.clone(),
            columns,
        };
        edited.parse_noodles()?;
        *self = edited;
        Ok(())
    }

    pub(super) fn apply_fai(&mut self, path: &Path) -> Result<()> {
        let fai = Fai::read(path)?;
        let edited = Self {
            metadata: rewrite_contigs(&self.metadata, &fai)?,
            columns: self.columns.clone(),
        };
        edited.parse_noodles()?;
        *self = edited;
        Ok(())
    }
}

fn parse_columns(line: &str) -> Result<Vec<String>> {
    let columns: Vec<_> = line.split('\t').map(str::to_owned).collect();
    const FIXED: [&str; 8] = [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
    ];
    if columns.len() < FIXED.len() || columns[..FIXED.len()] != FIXED {
        return Err(invalid("#CHROM line has invalid fixed columns"));
    }
    if columns.len() == 9 || (columns.len() > 9 && columns[8] != "FORMAT") {
        return Err(invalid(
            "#CHROM sample columns require FORMAT and at least one sample",
        ));
    }
    if columns.len() > 9 {
        let mut names = HashSet::with_capacity(columns.len() - 9);
        for name in &columns[9..] {
            if name.is_empty() {
                return Err(invalid("VCF sample names cannot be empty"));
            }
            if !names.insert(name) {
                return Err(invalid(format!("duplicate VCF sample name: {name}")));
            }
        }
    }
    Ok(columns)
}

fn validate_line(line: &str) -> Result<()> {
    if let Some(character) = line
        .chars()
        .find(|character| character.is_control() && *character != '\t')
    {
        return Err(invalid(format!(
            "VCF header contains disallowed control character U+{:04X}",
            u32::from(character)
        )));
    }
    Ok(())
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use crate::reheader::samples::SampleSource;

    use super::HeaderText;

    const HEADER_TWO_SAMPLES: &str = "##fileformat=VCFv4.3\n\
##source=original\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

    const REPLACEMENT_CRLF: &str = "##fileformat=VCFv4.3\r\n\
##source=changed\r\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tR1\tR2\r\n";

    fn fixture(source: &[u8]) -> tempfile::NamedTempFile {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(source).unwrap();
        file
    }

    #[test]
    fn replacement_preserves_unknown_lines_and_normalizes_newlines() {
        let old = HeaderText::parse(HEADER_TWO_SAMPLES.as_bytes()).unwrap();
        let input = fixture(REPLACEMENT_CRLF.as_bytes());
        let replacement = HeaderText::replace_from(input.path(), old.sample_names().len()).unwrap();
        assert_eq!(
            replacement.render(),
            b"##fileformat=VCFv4.3\n##source=changed\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tR1\tR2\n"
        );
    }

    #[test]
    fn replacement_rejects_a_different_sample_count() {
        let input = fixture(
            b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tR1\n",
        );
        let error = HeaderText::replace_from(input.path(), 2).unwrap_err();
        assert!(error.to_string().contains("2 samples"), "{error}");
    }

    #[test]
    fn header_shape_is_strict() {
        for source in [
            "",
            "##source=x\n##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "##fileformat=VCFv4.3\n##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n##source=late\n",
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tWRONG\n",
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tGT\tS1\n",
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\n",
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS1\n",
        ] {
            assert!(HeaderText::parse(source.as_bytes()).is_err(), "{source:?}");
        }
    }

    #[test]
    fn rejects_non_utf8_and_disallowed_controls() {
        let mut non_utf8 = HEADER_TWO_SAMPLES.as_bytes().to_vec();
        non_utf8[2] = 0xff;
        assert!(HeaderText::parse(&non_utf8).is_err());

        for control in [0, 0x0b, 0x0c, 0x7f] {
            let mut source = b"##fileformat=VCFv4.3\n##source=x\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec();
            source[31] = control;
            assert!(HeaderText::parse(&source).is_err(), "control={control}");
        }
    }

    #[test]
    fn exposes_samples_contigs_and_a_noodles_header() {
        let header = HeaderText::parse(HEADER_TWO_SAMPLES.as_bytes()).unwrap();
        assert_eq!(header.sample_names(), ["S1", "S2"]);
        assert_eq!(header.contig_count(), 1);
        let parsed = header.parse_noodles().unwrap();
        assert_eq!(parsed.sample_names().len(), 2);
        assert_eq!(parsed.contigs().len(), 1);
    }

    #[test]
    fn sample_edits_replace_only_the_terminal_columns() {
        let mut header = HeaderText::parse(HEADER_TWO_SAMPLES.as_bytes()).unwrap();
        header
            .apply_samples(&SampleSource::List("Tumor,Normal".to_owned()))
            .unwrap();
        assert_eq!(
            header.render(),
            b"##fileformat=VCFv4.3\n##source=original\n##contig=<ID=chr1,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTumor\tNormal\n"
        );
    }
}
