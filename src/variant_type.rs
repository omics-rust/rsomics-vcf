const REF: u32 = 0;
const SNP: u32 = 1 << 0;
const MNP: u32 = 1 << 1;
const INDEL: u32 = 1 << 2;
const OTHER: u32 = 1 << 3;
const BND: u32 = 1 << 4;
const OVERLAP: u32 = 1 << 5;

pub(crate) fn write(output: &mut Vec<u8>, reference: &[u8], alternates: &[u8]) {
    let mut mask = REF;
    for alternate in alternates.split(|value| *value == b',') {
        mask |= classify(reference, alternate);
    }

    let mut fields = 0;
    for (bit, name) in [
        (SNP, b"SNP".as_slice()),
        (MNP, b"MNP".as_slice()),
        (INDEL, b"INDEL".as_slice()),
        (OTHER, b"OTHER".as_slice()),
        (BND, b"BND".as_slice()),
        (OVERLAP, b"OVERLAP".as_slice()),
    ] {
        if mask & bit != 0 {
            if fields > 0 {
                output.push(b',');
            }
            output.extend_from_slice(name);
            fields += 1;
        }
    }
    if fields == 0 {
        output.extend_from_slice(b"REF");
    }
}

fn classify(reference: &[u8], alternate: &[u8]) -> u32 {
    if alternate == b"*" {
        return OVERLAP;
    }
    if reference.len() == 1 && alternate.len() == 1 {
        return if alternate[0] == b'.'
            || reference[0].eq_ignore_ascii_case(&alternate[0])
            || alternate[0] == b'X'
        {
            REF
        } else {
            SNP
        };
    }
    if alternate.first() == Some(&b'<') {
        return if alternate.starts_with(b"<X>")
            || alternate.starts_with(b"<*>")
            || alternate == b"<NON_REF>"
        {
            REF
        } else {
            OTHER
        };
    }
    if matches!(alternate.first(), Some(b']') | Some(b'[')) {
        return BND;
    }

    let (mut reference_start, mut alternate_start) = (0, 0);
    while reference_start < reference.len()
        && alternate_start < alternate.len()
        && reference[reference_start].eq_ignore_ascii_case(&alternate[alternate_start])
    {
        reference_start += 1;
        alternate_start += 1;
    }
    if alternate_start < alternate.len() && reference_start == reference.len() {
        return if matches!(alternate.last(), Some(b']') | Some(b'[')) {
            BND
        } else {
            INDEL
        };
    }
    if reference_start < reference.len() && alternate_start == alternate.len() {
        return INDEL;
    }
    if reference_start == reference.len() && alternate_start == alternate.len() {
        return REF;
    }

    let (mut reference_end, mut alternate_end) = (reference.len() - 1, alternate.len() - 1);
    if matches!(alternate[alternate_end], b']' | b'[') {
        return BND;
    }
    while reference_end > reference_start
        && alternate_end > alternate_start
        && reference[reference_end].eq_ignore_ascii_case(&alternate[alternate_end])
    {
        reference_end -= 1;
        alternate_end -= 1;
    }
    if alternate_end == alternate_start {
        return if reference_end == reference_start {
            SNP
        } else if reference[reference_end].eq_ignore_ascii_case(&alternate[alternate_end]) {
            INDEL
        } else {
            OTHER
        };
    }
    if reference_end == reference_start {
        return if reference[reference_end].eq_ignore_ascii_case(&alternate[alternate_end]) {
            INDEL
        } else {
            OTHER
        };
    }
    if reference_end - reference_start == alternate_end - alternate_start {
        MNP
    } else {
        OTHER
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kind(reference: &str, alternates: &str) -> String {
        let mut output = Vec::new();
        write(&mut output, reference.as_bytes(), alternates.as_bytes());
        String::from_utf8(output).unwrap()
    }

    #[test]
    fn reports_current_bcftools_categories() {
        assert_eq!(kind("A", "."), "REF");
        assert_eq!(kind("A", "G"), "SNP");
        assert_eq!(kind("AC", "GT"), "MNP");
        assert_eq!(kind("A", "AT"), "INDEL");
        assert_eq!(kind("A", "<DEL>"), "OTHER");
        assert_eq!(kind("A", "A]chr2:4]"), "BND");
        assert_eq!(kind("A", "*"), "OVERLAP");
        assert_eq!(kind("A", "C,GG"), "SNP,OTHER");
    }
}
