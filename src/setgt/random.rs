#[derive(Clone, Copy, Debug)]
pub(super) struct Random48 {
    state: u64,
}

impl Random48 {
    pub fn new(seed: i64) -> Self {
        Self {
            state: (u64::from(seed as u32) << 16) | 0x330e,
        }
    }

    pub fn next(&mut self) -> f64 {
        self.state = self.state.wrapping_mul(0x5deece66d).wrapping_add(0xb) & ((1_u64 << 48) - 1);
        self.state as f64 / (1_u64 << 48) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seed_seven_matches_htslib_1_24_drand48_bits() {
        let mut random = Random48::new(7);
        for expected in [
            0x3fd10d6bf5d44040,
            0x3fe5d33b8c0c6f00,
            0x3fd0fdcc420a88c0,
            0x3fc086b44cb17900,
        ] {
            assert_eq!(random.next().to_bits(), expected);
        }
    }

    #[test]
    fn signed_seeds_use_the_low_32_twos_complement_bits() {
        let mut signed = Random48::new(-1);
        let mut low_bits = Random48::new(u32::MAX as i64);
        for _ in 0..4 {
            assert_eq!(signed.next().to_bits(), low_bits.next().to_bits());
        }
    }
}
