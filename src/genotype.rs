mod counts;
mod edit;

pub(crate) use counts::{InfoPolicy, allele_counts, reconcile_ac_an};
pub(crate) use edit::{Change, MissingPolicy, edit_selected};

#[cfg(test)]
pub(crate) use counts::AlleleCounts;

#[cfg(test)]
mod tests {
    use noodles_vcf::{
        self as vcf,
        header::record::value::{
            Map,
            map::{
                Info,
                info::{Number, Type},
            },
        },
        variant::{
            RecordBuf,
            record::samples::series::value::genotype::Phasing,
            record_buf::{
                Samples,
                info::field::{Value as InfoValue, value::Array},
                samples::sample::{
                    Value as SampleValue,
                    value::{Genotype, genotype::Allele},
                },
            },
        },
    };

    use super::*;

    fn fixture(samples: &[&str], info: &str) -> (vcf::Header, RecordBuf) {
        let names = (1..=samples.len())
            .map(|index| format!("S{index}"))
            .collect::<Vec<_>>()
            .join("\t");
        let header: vcf::Header = format!(
            "##fileformat=VCFv4.3\n\
             ##INFO=<ID=AC,Number=A,Type=Integer,Description=\"allele count\">\n\
             ##INFO=<ID=AN,Number=1,Type=Integer,Description=\"allele number\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{names}\n"
        )
        .parse()
        .unwrap();
        let line = format!(
            "chr1\t1\t.\tA\tC,G\t10\tPASS\t{info}\tGT\t{}",
            samples.join("\t")
        );
        let raw = vcf::Record::try_from(line.as_bytes()).unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        (header, record)
    }

    fn no_genotype_fixture() -> RecordBuf {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            .parse()
            .unwrap();
        let raw =
            vcf::Record::try_from(b"chr1\t1\t.\tA\tC\t10\tPASS\t.\tDP\t8".as_slice()).unwrap();
        RecordBuf::try_from_variant_record(&header, &raw).unwrap()
    }

    fn alleles(record: &RecordBuf, sample: usize) -> Vec<(Option<usize>, Phasing)> {
        let sample = record.samples().get_index(sample).unwrap();
        let Some(Some(SampleValue::Genotype(genotype))) = sample.get("GT") else {
            return Vec::new();
        };
        genotype
            .as_ref()
            .iter()
            .map(|allele| (allele.position(), allele.phasing()))
            .collect()
    }

    fn integer_info(record: &RecordBuf, name: &str) -> Option<i32> {
        match record.info().get(name) {
            Some(Some(InfoValue::Integer(value))) => Some(*value),
            _ => None,
        }
    }

    fn integer_array_info(record: &RecordBuf, name: &str) -> Option<Vec<Option<i32>>> {
        match record.info().get(name) {
            Some(Some(InfoValue::Array(Array::Integer(values)))) => Some(values.clone()),
            _ => None,
        }
    }

    #[test]
    fn selected_edits_support_arbitrary_ploidy_and_count_allele_changes() {
        let (_, mut record) = fixture(&["0|1", "1/1", "0/1/2"], "AC=3,1;AN=7");
        let change = edit_selected(
            &mut record,
            &[false, true, true],
            MissingPolicy::Error,
            |_, genotype| {
                Ok(genotype
                    .as_ref()
                    .iter()
                    .map(|_| Allele::new(Some(0), Phasing::Unphased))
                    .collect())
            },
        )
        .unwrap();

        assert_eq!(
            change,
            Change {
                genotypes: 2,
                alleles: 4,
            }
        );
        assert_eq!(
            alleles(&record, 0),
            [(Some(0), Phasing::Phased), (Some(1), Phasing::Phased)]
        );
        assert_eq!(
            alleles(&record, 1),
            [(Some(0), Phasing::Unphased), (Some(0), Phasing::Unphased)]
        );
        assert_eq!(alleles(&record, 2).len(), 3);
    }

    #[test]
    fn ploidy_growth_shrink_and_phase_changes_are_counted() {
        let (_, mut record) = fixture(&["0/1", "0/1/2"], "AC=2,1;AN=5");
        let change = edit_selected(
            &mut record,
            &[true, true],
            MissingPolicy::Error,
            |sample, _| {
                let genotype: Genotype = if sample == 0 {
                    "0/1|2".parse().unwrap()
                } else {
                    [Allele::new(Some(0), Phasing::Unphased)]
                        .into_iter()
                        .collect()
                };
                Ok(genotype)
            },
        )
        .unwrap();

        assert_eq!(
            change,
            Change {
                genotypes: 2,
                alleles: 3,
            }
        );
        assert_eq!(alleles(&record, 0).len(), 3);
        assert_eq!(alleles(&record, 1), [(Some(0), Phasing::Unphased)]);
    }

    #[test]
    fn phase_only_edits_count_each_changed_allele() {
        let (_, mut record) = fixture(&["0|1"], "AC=1,0;AN=2");
        let change = edit_selected(&mut record, &[true], MissingPolicy::Error, |_, _| {
            Ok("0/1".parse().unwrap())
        })
        .unwrap();
        assert_eq!(
            change,
            Change {
                genotypes: 1,
                alleles: 2,
            }
        );
    }

    #[test]
    fn edit_validation_is_transactional() {
        let (_, mut record) = fixture(&["0/1", "1/1"], "AC=3,0;AN=4");
        let original = record.clone();
        assert!(
            edit_selected(
                &mut record,
                &[true],
                MissingPolicy::Error,
                |_, genotype| Ok(genotype.clone())
            )
            .is_err()
        );
        assert_eq!(record, original);

        let original = record.clone();
        assert!(
            edit_selected(
                &mut record,
                &[true, true],
                MissingPolicy::Error,
                |sample, _| {
                    if sample == 1 {
                        Err(rsomics_common::RsomicsError::InvalidInput(
                            "injected edit failure".to_owned(),
                        ))
                    } else {
                        Ok("./.".parse().unwrap())
                    }
                }
            )
            .is_err()
        );
        assert_eq!(record, original);

        let mut no_genotype = no_genotype_fixture();
        assert!(
            edit_selected(
                &mut no_genotype,
                &[true],
                MissingPolicy::Error,
                |_, genotype| Ok(genotype.clone())
            )
            .is_err()
        );
        assert_eq!(
            edit_selected(
                &mut no_genotype,
                &[true],
                MissingPolicy::Ignore,
                |_, genotype| Ok(genotype.clone())
            )
            .unwrap(),
            Change::default()
        );
    }

    #[test]
    fn missing_values_follow_policy_and_wrong_types_fail() {
        let (_, mut record) = fixture(&[".", "0/1"], "AC=1,0;AN=2");
        let mut visited = Vec::new();
        edit_selected(
            &mut record,
            &[true, true],
            MissingPolicy::Ignore,
            |sample, genotype| {
                visited.push(sample);
                Ok(genotype.clone())
            },
        )
        .unwrap();
        assert_eq!(visited, [1]);
        assert!(
            edit_selected(
                &mut record,
                &[true, false],
                MissingPolicy::Error,
                |_, genotype| Ok(genotype.clone())
            )
            .is_err()
        );

        let keys = record.samples().keys().clone();
        let genotype_index = keys.as_ref().get_index_of("GT").unwrap();
        let mut rows = record
            .samples()
            .values()
            .map(|sample| sample.values().to_vec())
            .collect::<Vec<_>>();
        rows[1][genotype_index] = Some(SampleValue::Integer(1));
        *record.samples_mut() = Samples::new(keys, rows);
        assert!(
            edit_selected(
                &mut record,
                &[false, true],
                MissingPolicy::Ignore,
                |_, genotype| Ok(genotype.clone())
            )
            .is_err()
        );
    }

    #[test]
    fn allele_counts_include_reference_and_reject_out_of_range_positions() {
        let (_, record) = fixture(&["0/1", "2/.", "0/0/0"], "AC=1,1;AN=6");
        assert_eq!(
            allele_counts(&record).unwrap(),
            AlleleCounts {
                counts: vec![4, 1, 1],
                total: 6,
            }
        );

        let (_, record) = fixture(&["0/3"], "AC=0,0;AN=2");
        assert!(allele_counts(&record).is_err());
    }

    #[test]
    fn reconciliation_updates_present_valid_tags_and_keeps_absent_tags_absent() {
        let (header, mut record) = fixture(&["0/0", "0/0/0"], "AC=8,9;AN=1");
        reconcile_ac_an(&header, &mut record, InfoPolicy::Strict).unwrap();
        assert_eq!(
            integer_array_info(&record, "AC"),
            Some(vec![Some(0), Some(0)])
        );
        assert_eq!(integer_info(&record, "AN"), Some(5));

        record.info_mut().as_mut().shift_remove("AN");
        reconcile_ac_an(&header, &mut record, InfoPolicy::Strict).unwrap();
        assert!(!record.info().as_ref().contains_key("AN"));
    }

    #[test]
    fn strict_reconciliation_rejects_invalid_definitions_and_values() {
        let (mut header, mut record) = fixture(&["0/1", "1/1"], "AC=9,9;AN=9");
        header.infos_mut().insert(
            "AC".into(),
            Map::<Info>::new(Number::Count(1), Type::Integer, "wrong number"),
        );
        let original = record.info().clone();
        assert!(reconcile_ac_an(&header, &mut record, InfoPolicy::Strict).is_err());
        assert_eq!(record.info(), &original);

        let (header, mut record) = fixture(&["0/1", "1/1"], "AC=9,9;AN=9");
        record
            .info_mut()
            .insert("AC".into(), Some(InfoValue::Integer(9)));
        let original = record.info().clone();
        assert!(reconcile_ac_an(&header, &mut record, InfoPolicy::Strict).is_err());
        assert_eq!(record.info(), &original);
    }

    #[test]
    fn best_effort_reconciliation_preserves_each_invalid_tag_independently() {
        let (mut header, mut record) = fixture(&["0/1", "1/1"], "AC=9,9;AN=9");
        header.infos_mut().insert(
            "AC".into(),
            Map::<Info>::new(Number::Count(1), Type::Integer, "wrong number"),
        );
        reconcile_ac_an(&header, &mut record, InfoPolicy::BestEffort).unwrap();
        assert_eq!(
            integer_array_info(&record, "AC"),
            Some(vec![Some(9), Some(9)])
        );
        assert_eq!(integer_info(&record, "AN"), Some(4));

        let (header, mut record) = fixture(&["0/1", "1/1"], "AC=9,9;AN=9");
        record
            .info_mut()
            .insert("AC".into(), Some(InfoValue::Integer(9)));
        reconcile_ac_an(&header, &mut record, InfoPolicy::BestEffort).unwrap();
        assert_eq!(record.info().get("AC"), Some(Some(&InfoValue::Integer(9))));
        assert_eq!(integer_info(&record, "AN"), Some(4));
    }
}
