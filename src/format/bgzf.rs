use std::io::{self, Cursor, Read, Write};

pub(crate) const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

#[derive(Debug, Eq, PartialEq)]
pub(crate) enum Frame {
    Data(Vec<u8>),
    Eof,
}

pub(crate) struct FrameReader<R> {
    inner: R,
    finished: bool,
}

impl<R> FrameReader<R> {
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            finished: false,
        }
    }

    pub(crate) fn into_inner(self) -> R {
        self.inner
    }
}

impl<R: Read> FrameReader<R> {
    pub(crate) fn next(&mut self) -> io::Result<Option<Frame>> {
        if self.finished {
            return Ok(None);
        }

        let mut header = [0; 18];
        let Some(first) = read_byte(&mut self.inner)? else {
            return Ok(None);
        };
        header[0] = first;
        read_complete(&mut self.inner, &mut header[1..], "partial BGZF header")?;
        validate_header(&header)?;

        let length = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;
        if !(EOF_BLOCK.len()..=65_536).contains(&length) {
            return Err(invalid(format!("invalid BGZF block size: {length}")));
        }
        let mut raw = Vec::with_capacity(length);
        raw.extend_from_slice(&header);
        raw.resize(length, 0);
        read_complete(
            &mut self.inner,
            &mut raw[header.len()..],
            "partial BGZF payload",
        )?;

        if raw == EOF_BLOCK {
            if read_byte(&mut self.inner)?.is_some() {
                return Err(invalid("bytes follow the canonical BGZF EOF block"));
            }
            self.finished = true;
            Ok(Some(Frame::Eof))
        } else {
            Ok(Some(Frame::Data(raw)))
        }
    }

    pub(crate) fn copy_through_eof<W: Write>(&mut self, mut output: W) -> io::Result<u64> {
        let mut copied = 0u64;
        loop {
            let raw = match self.next()? {
                Some(Frame::Data(raw)) => raw,
                Some(Frame::Eof) => EOF_BLOCK.to_vec(),
                None => return Err(invalid("canonical BGZF EOF block is missing")),
            };
            output.write_all(&raw)?;
            copied = copied
                .checked_add(raw.len() as u64)
                .ok_or_else(|| invalid("copied BGZF length exceeds u64"))?;
            if raw == EOF_BLOCK {
                return Ok(copied);
            }
        }
    }
}

pub(crate) struct ValidatedReader<R> {
    frames: FrameReader<R>,
    current: Cursor<Vec<u8>>,
    complete: bool,
}

impl<R> ValidatedReader<R> {
    pub(crate) fn new(inner: R) -> Self {
        Self {
            frames: FrameReader::new(inner),
            current: Cursor::new(Vec::new()),
            complete: false,
        }
    }
}

impl<R: Read> Read for ValidatedReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        if buffer.is_empty() {
            return Ok(0);
        }
        loop {
            let read = self.current.read(buffer)?;
            if read != 0 {
                return Ok(read);
            }
            if self.complete {
                return Ok(0);
            }
            let raw = match self.frames.next()? {
                Some(Frame::Data(raw)) => raw,
                Some(Frame::Eof) => {
                    self.complete = true;
                    EOF_BLOCK.to_vec()
                }
                None => return Err(invalid("canonical BGZF EOF block is missing")),
            };
            self.current = Cursor::new(raw);
        }
    }
}

fn validate_header(header: &[u8; 18]) -> io::Result<()> {
    if header[..4] != [0x1f, 0x8b, 0x08, 0x04] {
        return Err(invalid("invalid BGZF gzip header"));
    }
    if header[10..16] != [0x06, 0x00, 0x42, 0x43, 0x02, 0x00] {
        return Err(invalid("invalid BGZF extra field"));
    }
    Ok(())
}

fn read_byte(reader: &mut impl Read) -> io::Result<Option<u8>> {
    let mut byte = [0];
    loop {
        match reader.read(&mut byte) {
            Ok(0) => return Ok(None),
            Ok(_) => return Ok(Some(byte[0])),
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
}

fn read_complete(reader: &mut impl Read, buffer: &mut [u8], message: &str) -> io::Result<()> {
    let mut offset = 0;
    while offset < buffer.len() {
        match reader.read(&mut buffer[offset..]) {
            Ok(0) => return Err(invalid(message)),
            Ok(length) => offset += length,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(())
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

#[cfg(test)]
mod tests {
    use std::io::{self, Write};

    use super::{EOF_BLOCK, Frame, FrameReader};

    fn bgzf_fixture(source: &[u8]) -> Vec<u8> {
        let mut writer = noodles_bgzf::io::Writer::new(Vec::new());
        writer.write_all(source).unwrap();
        writer.finish().unwrap()
    }

    #[test]
    fn copies_complete_frames_through_one_canonical_eof() {
        let input = bgzf_fixture(b"body");
        let mut output = Vec::new();
        let copied = FrameReader::new(input.as_slice())
            .copy_through_eof(&mut output)
            .unwrap();
        assert_eq!(copied, input.len() as u64);
        assert_eq!(output, input);
    }

    #[test]
    fn exposes_data_and_eof_frames_in_order() {
        let input = bgzf_fixture(b"body");
        let mut reader = FrameReader::new(input.as_slice());
        assert!(matches!(reader.next().unwrap(), Some(Frame::Data(_))));
        assert!(matches!(reader.next().unwrap(), Some(Frame::Eof)));
        assert!(reader.next().unwrap().is_none());
    }

    #[test]
    fn rejects_bytes_or_a_second_frame_after_eof() {
        for suffix in [b"trailing".as_slice(), EOF_BLOCK.as_slice()] {
            let mut input = bgzf_fixture(b"body");
            input.extend_from_slice(suffix);
            assert_eq!(
                FrameReader::new(input.as_slice())
                    .copy_through_eof(io::sink())
                    .unwrap_err()
                    .kind(),
                io::ErrorKind::InvalidData
            );
        }
    }

    #[test]
    fn rejects_partial_headers_and_payloads() {
        let input = bgzf_fixture(b"body");
        for end in [1, 17, input.len() - EOF_BLOCK.len() - 1] {
            assert_eq!(
                FrameReader::new(&input[..end])
                    .copy_through_eof(io::sink())
                    .unwrap_err()
                    .kind(),
                io::ErrorKind::InvalidData,
                "end={end}"
            );
        }
    }

    #[test]
    fn rejects_invalid_frame_headers_and_sizes() {
        let input = bgzf_fixture(b"body");
        for (index, value) in [(0, 0), (3, 0), (10, 0), (12, b'X'), (14, 0)] {
            let mut invalid = input.clone();
            invalid[index] = value;
            assert_eq!(
                FrameReader::new(invalid.as_slice())
                    .copy_through_eof(io::sink())
                    .unwrap_err()
                    .kind(),
                io::ErrorKind::InvalidData,
                "index={index}"
            );
        }

        let mut invalid = input;
        invalid[16..18].copy_from_slice(&0u16.to_le_bytes());
        assert_eq!(
            FrameReader::new(invalid.as_slice())
                .copy_through_eof(io::sink())
                .unwrap_err()
                .kind(),
            io::ErrorKind::InvalidData
        );
    }

    #[test]
    fn rejects_a_stream_without_the_canonical_eof() {
        let mut input = bgzf_fixture(b"body");
        input.truncate(input.len() - EOF_BLOCK.len());
        assert_eq!(
            FrameReader::new(input.as_slice())
                .copy_through_eof(io::sink())
                .unwrap_err()
                .kind(),
            io::ErrorKind::InvalidData
        );
    }

    #[test]
    fn accepts_a_noncanonical_empty_data_frame_before_eof() {
        let mut empty = EOF_BLOCK;
        empty[4] = 1;
        let mut input = empty.to_vec();
        input.extend_from_slice(&EOF_BLOCK);
        let mut reader = FrameReader::new(input.as_slice());
        assert!(matches!(reader.next().unwrap(), Some(Frame::Data(frame)) if frame == empty));
        assert!(matches!(reader.next().unwrap(), Some(Frame::Eof)));
    }
}
