pub(super) fn infer_ploidy(alleles: usize, values: usize) -> Option<usize> {
    (1..=64).find(|&ploidy| combinations(alleles + ploidy - 1, ploidy) == Some(values))
}

pub(super) fn combinations(n: usize, k: usize) -> Option<usize> {
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
