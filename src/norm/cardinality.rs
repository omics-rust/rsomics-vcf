pub(crate) fn infer_ploidy(alleles: usize, values: usize) -> Option<usize> {
    (1..=64).find(|&ploidy| {
        alleles
            .checked_add(ploidy - 1)
            .and_then(|width| combinations(width, ploidy))
            == Some(values)
    })
}

pub(crate) fn combinations(n: usize, k: usize) -> Option<usize> {
    if k > n {
        return Some(0);
    }
    let k = k.min(n - k);
    (1..=k).try_fold(1usize, |value, divisor| {
        value
            .checked_mul(n - k + divisor)
            .map(|product| product / divisor)
    })
}

pub(crate) fn genotype_index(genotype: &[usize]) -> Option<usize> {
    genotype
        .iter()
        .enumerate()
        .try_fold(0usize, |index, (i, allele)| {
            allele
                .checked_add(i)
                .and_then(|width| combinations(width, i + 1))
                .and_then(|offset| index.checked_add(offset))
        })
}

pub(crate) fn visit_genotypes(alleles: usize, ploidy: usize, mut visit: impl FnMut(&[usize])) {
    visit_genotypes_inner(alleles, ploidy, 0, &mut Vec::new(), &mut visit);
}

fn visit_genotypes_inner(
    alleles: usize,
    ploidy: usize,
    minimum: usize,
    genotype: &mut Vec<usize>,
    visit: &mut impl FnMut(&[usize]),
) {
    if genotype.len() == ploidy {
        visit(genotype);
        return;
    }
    for allele in minimum..alleles {
        genotype.push(allele);
        visit_genotypes_inner(alleles, ploidy, allele, genotype, visit);
        genotype.pop();
    }
}
