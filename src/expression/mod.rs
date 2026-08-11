mod bind;
mod evaluate;
mod syntax;
mod value;

use std::fmt;

use noodles_vcf::{self as vcf, variant::RecordBuf};

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Compiled(bind::BoundExpression);

impl Compiled {
    pub(crate) fn bind(source: &str, header: &vcf::Header) -> Result<Self, CompileError> {
        let expression = syntax::parse(source).map_err(CompileError::Parse)?;
        bind::bind(expression, header)
            .map(Self)
            .map_err(CompileError::Bind)
    }

    pub(crate) fn evaluate(
        &self,
        header: &vcf::Header,
        record: &RecordBuf,
    ) -> Result<evaluate::Truth, evaluate::EvaluateError> {
        evaluate::evaluate(&self.0, header, record)
    }
}

#[derive(Debug)]
pub(crate) enum CompileError {
    Parse(syntax::ParseError),
    Bind(bind::BindError),
}

impl fmt::Display for CompileError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(error) => error.fmt(formatter),
            Self::Bind(error) => error.fmt(formatter),
        }
    }
}
