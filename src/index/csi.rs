use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{
    self as csi, BinningIndex,
    binning_index::{
        self,
        index::{
            Header, ReferenceSequence,
            reference_sequence::{Metadata, index::BinnedIndex},
        },
    },
};

const MAX_DEPTH: u8 = 9;

pub(super) fn settings(mut min_shift: u8, mut depth: u8, max_length: Option<usize>) -> (u8, u8) {
    let Some(max_length) = max_length else {
        depth = if min_shift < 10 {
            MAX_DEPTH
        } else if min_shift < 25 {
            MAX_DEPTH - (min_shift - 10) / 3
        } else {
            4
        };
        return (min_shift, depth);
    };
    let required = max_length as u128 + 256;
    while depth < MAX_DEPTH && required > capacity(min_shift, depth) {
        depth += 1;
    }
    while required > capacity(min_shift, depth) {
        min_shift += 1;
    }
    (min_shift, depth)
}

fn capacity(min_shift: u8, depth: u8) -> u128 {
    1u128 << (min_shift as u32 + 3 * depth as u32)
}

pub(super) struct LinearIndex {
    min_shift: u8,
    depth: u8,
    offsets: Vec<Option<bgzf::VirtualPosition>>,
}

impl LinearIndex {
    pub(super) fn new(min_shift: u8, depth: u8) -> Self {
        Self {
            min_shift,
            depth,
            offsets: Vec::new(),
        }
    }

    pub(super) fn insert(&mut self, start: Position, end: Position, offset: bgzf::VirtualPosition) {
        let begin = (usize::from(start) - 1) >> self.min_shift;
        let end = (usize::from(end) - 1) >> self.min_shift;
        if self.offsets.len() <= end {
            self.offsets.resize(end + 1, None);
        }
        for target in &mut self.offsets[begin..=end] {
            target.get_or_insert(offset);
        }
    }

    fn finish(&mut self) {
        for i in (0..self.offsets.len().saturating_sub(1)).rev() {
            if self.offsets[i].is_none() {
                self.offsets[i] = self.offsets[i + 1];
            }
        }
    }

    fn loffset(&self, bin: usize) -> bgzf::VirtualPosition {
        let bottom = bin_bottom(bin, self.depth);
        self.offsets
            .get(bottom)
            .copied()
            .flatten()
            .unwrap_or_default()
    }
}

fn bin_bottom(bin: usize, depth: u8) -> usize {
    let level = bin_level(bin);
    (bin - bin_first(level)) << ((u32::from(depth) - level) * 3)
}

fn bin_level(mut bin: usize) -> u32 {
    let mut level = 0;
    while bin != 0 {
        bin = (bin - 1) >> 3;
        level += 1;
    }
    level
}

fn bin_first(level: u32) -> usize {
    ((1usize << (3 * level)) - 1) / 7
}

pub(super) fn apply_loffsets(
    index: csi::Index,
    linear: &mut [LinearIndex],
    header: Option<Header>,
) -> csi::Index {
    let reference_sequences = index
        .reference_sequences()
        .iter()
        .zip(linear)
        .map(|(reference_sequence, linear)| {
            linear.finish();
            let loffsets = reference_sequence
                .bins()
                .keys()
                .map(|&bin| (bin, linear.loffset(bin)))
                .collect::<BinnedIndex>();
            ReferenceSequence::new(
                reference_sequence.bins().clone(),
                loffsets,
                metadata(reference_sequence).cloned(),
            )
        })
        .collect();
    let mut builder = csi::Index::builder()
        .set_min_shift(index.min_shift())
        .set_depth(index.depth())
        .set_reference_sequences(reference_sequences);
    if let Some(header) = header.or_else(|| index.header().cloned()) {
        builder = builder.set_header(header);
    }
    if let Some(count) = index.unplaced_unmapped_record_count() {
        builder = builder.set_unplaced_unmapped_record_count(count);
    }
    builder.build()
}

fn metadata<I>(reference: &ReferenceSequence<I>) -> Option<&Metadata>
where
    I: Default,
    ReferenceSequence<I>: binning_index::ReferenceSequence,
{
    use binning_index::ReferenceSequence as _;
    reference.metadata()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn settings_match_htslib_defaults() {
        assert_eq!(settings(14, 6, Some(248_956_422)), (14, 6));
        assert_eq!(settings(14, 0, Some(248_956_422)), (14, 5));
        assert_eq!(settings(14, 6, None), (14, 8));
    }
}
