use std::fmt;
use std::path::PathBuf;

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Expression {
    Value(Value),
    Unary {
        operator: UnaryOperator,
        expression: Box<Self>,
    },
    Binary {
        operator: BinaryOperator,
        left: Box<Self>,
        right: Box<Self>,
    },
    Function {
        name: String,
        arguments: Vec<Self>,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Value {
    Number(f64),
    String(String),
    Missing,
    File(PathBuf),
    Field(Field),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum UnaryOperator {
    Negate,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum BinaryOperator {
    Add,
    Subtract,
    Multiply,
    Divide,
    Modulo,
    Equal,
    NotEqual,
    Less,
    LessEqual,
    Greater,
    GreaterEqual,
    Regex,
    NotRegex,
    SampleAnd,
    SiteAnd,
    SampleOr,
    SiteOr,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Namespace {
    Info,
    Format,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Field {
    pub namespace: Option<Namespace>,
    pub name: String,
    pub subscript: Option<Subscript>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Subscript {
    Values(Selector),
    SampleValues { samples: Selector, values: Selector },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Selector {
    Any,
    Genotype,
    File(PathBuf),
    Indices(Vec<IndexSelector>),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct IndexSelector {
    pub start: usize,
    pub end: Option<usize>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ParseError {
    offset: usize,
    message: String,
}

impl ParseError {
    fn new(offset: usize, message: impl Into<String>) -> Self {
        Self {
            offset,
            message: message.into(),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "expression byte {}: {}",
            self.offset, self.message
        )
    }
}

impl std::error::Error for ParseError {}

#[derive(Clone, Debug, PartialEq)]
enum TokenKind {
    Number(f64),
    String(String),
    Identifier(String),
    File(PathBuf),
    LeftParen,
    RightParen,
    LeftBracket,
    RightBracket,
    Comma,
    Colon,
    Plus,
    Minus,
    Star,
    Slash,
    Percent,
    Equal,
    NotEqual,
    Less,
    LessEqual,
    Greater,
    GreaterEqual,
    Regex,
    NotRegex,
    Ampersand,
    DoubleAmpersand,
    Pipe,
    DoublePipe,
    End,
}

#[derive(Clone, Debug, PartialEq)]
struct Token {
    kind: TokenKind,
    offset: usize,
}

struct Lexer<'a> {
    source: &'a str,
    offset: usize,
}

impl<'a> Lexer<'a> {
    fn new(source: &'a str) -> Self {
        Self { source, offset: 0 }
    }

    fn next(&mut self) -> Result<Token, ParseError> {
        self.skip_whitespace();
        let offset = self.offset;
        let Some(byte) = self.peek() else {
            return Ok(Token {
                kind: TokenKind::End,
                offset,
            });
        };
        let kind = match byte {
            b'(' => self.one(TokenKind::LeftParen),
            b')' => self.one(TokenKind::RightParen),
            b'[' => self.one(TokenKind::LeftBracket),
            b']' => self.one(TokenKind::RightBracket),
            b',' => self.one(TokenKind::Comma),
            b':' => self.one(TokenKind::Colon),
            b'+' => self.one(TokenKind::Plus),
            b'-' => self.one(TokenKind::Minus),
            b'*' => self.one(TokenKind::Star),
            b'/' => self.one(TokenKind::Slash),
            b'%' => self.one(TokenKind::Percent),
            b'=' => {
                self.offset += 1;
                if self.peek() == Some(b'=') {
                    self.offset += 1;
                }
                TokenKind::Equal
            }
            b'!' => {
                self.offset += 1;
                match self.peek() {
                    Some(b'=') => {
                        self.offset += 1;
                        TokenKind::NotEqual
                    }
                    Some(b'~') => {
                        self.offset += 1;
                        TokenKind::NotRegex
                    }
                    _ => return Err(ParseError::new(offset, "expected = or ~ after !")),
                }
            }
            b'<' => {
                self.offset += 1;
                if self.peek() == Some(b'=') {
                    self.offset += 1;
                    TokenKind::LessEqual
                } else {
                    TokenKind::Less
                }
            }
            b'>' => {
                self.offset += 1;
                if self.peek() == Some(b'=') {
                    self.offset += 1;
                    TokenKind::GreaterEqual
                } else {
                    TokenKind::Greater
                }
            }
            b'~' => self.one(TokenKind::Regex),
            b'&' => {
                self.offset += 1;
                if self.peek() == Some(b'&') {
                    self.offset += 1;
                    TokenKind::DoubleAmpersand
                } else {
                    TokenKind::Ampersand
                }
            }
            b'|' => {
                self.offset += 1;
                if self.peek() == Some(b'|') {
                    self.offset += 1;
                    TokenKind::DoublePipe
                } else {
                    TokenKind::Pipe
                }
            }
            b'"' | b'\'' => self.string(byte)?,
            b'@' => self.file()?,
            b'.' if self.peek_next().is_some_and(|byte| byte.is_ascii_digit()) => self.number()?,
            b'0'..=b'9' => self.number()?,
            _ if !is_delimiter(byte) => self.identifier()?,
            _ => {
                return Err(ParseError::new(
                    offset,
                    format!("unexpected character {:?}", char::from(byte)),
                ));
            }
        };
        Ok(Token { kind, offset })
    }

    fn one(&mut self, kind: TokenKind) -> TokenKind {
        self.offset += 1;
        kind
    }

    fn peek(&self) -> Option<u8> {
        self.source.as_bytes().get(self.offset).copied()
    }

    fn peek_next(&self) -> Option<u8> {
        self.source.as_bytes().get(self.offset + 1).copied()
    }

    fn skip_whitespace(&mut self) {
        while self.peek().is_some_and(|byte| byte.is_ascii_whitespace()) {
            self.offset += 1;
        }
    }

    fn number(&mut self) -> Result<TokenKind, ParseError> {
        let start = self.offset;
        let mut decimal = false;
        while let Some(byte) = self.peek() {
            match byte {
                b'0'..=b'9' => self.offset += 1,
                b'.' if !decimal => {
                    decimal = true;
                    self.offset += 1;
                }
                _ => break,
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.offset += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.offset += 1;
            }
            let exponent = self.offset;
            while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                self.offset += 1;
            }
            if exponent == self.offset {
                return Err(ParseError::new(start, "number has no exponent digits"));
            }
        }
        let raw = &self.source[start..self.offset];
        let value = raw
            .parse()
            .map_err(|_| ParseError::new(start, "invalid number"))?;
        if !f64::is_finite(value) {
            return Err(ParseError::new(start, "number must be finite"));
        }
        Ok(TokenKind::Number(value))
    }

    fn string(&mut self, quote: u8) -> Result<TokenKind, ParseError> {
        let start = self.offset;
        self.offset += 1;
        let mut output = String::new();
        let mut segment = self.offset;
        loop {
            let Some(byte) = self.peek() else {
                return Err(ParseError::new(start, "unterminated string"));
            };
            if byte == quote {
                output.push_str(&self.source[segment..self.offset]);
                self.offset += 1;
                return Ok(TokenKind::String(output));
            }
            if byte == b'\\' {
                output.push_str(&self.source[segment..self.offset]);
                self.offset += 1;
                let Some(escaped) = self.peek() else {
                    return Err(ParseError::new(start, "unterminated string escape"));
                };
                output.push(char::from(escaped));
                self.offset += 1;
                segment = self.offset;
            } else {
                self.offset += 1;
            }
        }
    }

    fn file(&mut self) -> Result<TokenKind, ParseError> {
        let start = self.offset;
        self.offset += 1;
        let body = self.offset;
        while self.peek().is_some_and(|byte| !is_file_delimiter(byte)) {
            self.offset += 1;
        }
        if body == self.offset {
            return Err(ParseError::new(start, "file reference has no path"));
        }
        Ok(TokenKind::File(PathBuf::from(
            &self.source[body..self.offset],
        )))
    }

    fn identifier(&mut self) -> Result<TokenKind, ParseError> {
        let start = self.offset;
        while self.peek().is_some_and(|byte| !is_delimiter(byte)) {
            self.offset += 1;
        }
        if start == self.offset {
            return Err(ParseError::new(start, "empty identifier"));
        }
        Ok(TokenKind::Identifier(
            self.source[start..self.offset].to_owned(),
        ))
    }
}

fn is_delimiter(byte: u8) -> bool {
    byte.is_ascii_whitespace()
        || matches!(
            byte,
            b'(' | b')'
                | b'['
                | b']'
                | b','
                | b':'
                | b'+'
                | b'-'
                | b'*'
                | b'/'
                | b'%'
                | b'='
                | b'!'
                | b'<'
                | b'>'
                | b'~'
                | b'&'
                | b'|'
                | b'"'
                | b'\''
                | b'@'
        )
}

fn is_file_delimiter(byte: u8) -> bool {
    byte.is_ascii_whitespace()
        || matches!(
            byte,
            b'(' | b')' | b'[' | b']' | b',' | b':' | b'=' | b'!' | b'<' | b'>' | b'&' | b'|'
        )
}

struct Parser<'a> {
    lexer: Lexer<'a>,
    current: Token,
}

impl<'a> Parser<'a> {
    fn new(source: &'a str) -> Result<Self, ParseError> {
        let mut lexer = Lexer::new(source);
        let current = lexer.next()?;
        Ok(Self { lexer, current })
    }

    fn parse(mut self) -> Result<Expression, ParseError> {
        if self.current.kind == TokenKind::End {
            return Err(ParseError::new(0, "expression is empty"));
        }
        let expression = self.parse_site_or()?;
        if self.current.kind != TokenKind::End {
            return Err(ParseError::new(
                self.current.offset,
                "unexpected token after expression",
            ));
        }
        Ok(expression)
    }

    fn advance(&mut self) -> Result<Token, ParseError> {
        let next = self.lexer.next()?;
        Ok(std::mem::replace(&mut self.current, next))
    }

    fn parse_site_or(&mut self) -> Result<Expression, ParseError> {
        let mut expression = self.parse_sample_or()?;
        while self.current.kind == TokenKind::DoublePipe {
            self.advance()?;
            let right = self.parse_sample_or()?;
            expression = binary(BinaryOperator::SiteOr, expression, right);
        }
        Ok(expression)
    }

    fn parse_sample_or(&mut self) -> Result<Expression, ParseError> {
        let mut expression = self.parse_site_and()?;
        while self.current.kind == TokenKind::Pipe {
            self.advance()?;
            let right = self.parse_site_and()?;
            expression = binary(BinaryOperator::SampleOr, expression, right);
        }
        Ok(expression)
    }

    fn parse_site_and(&mut self) -> Result<Expression, ParseError> {
        let mut expression = self.parse_sample_and()?;
        while self.current.kind == TokenKind::DoubleAmpersand {
            self.advance()?;
            let right = self.parse_sample_and()?;
            expression = binary(BinaryOperator::SiteAnd, expression, right);
        }
        Ok(expression)
    }

    fn parse_sample_and(&mut self) -> Result<Expression, ParseError> {
        let mut expression = self.parse_comparison()?;
        while self.current.kind == TokenKind::Ampersand {
            self.advance()?;
            let right = self.parse_comparison()?;
            expression = binary(BinaryOperator::SampleAnd, expression, right);
        }
        Ok(expression)
    }

    fn parse_comparison(&mut self) -> Result<Expression, ParseError> {
        let left = self.parse_addition()?;
        let operator = match self.current.kind {
            TokenKind::Equal => BinaryOperator::Equal,
            TokenKind::NotEqual => BinaryOperator::NotEqual,
            TokenKind::Less => BinaryOperator::Less,
            TokenKind::LessEqual => BinaryOperator::LessEqual,
            TokenKind::Greater => BinaryOperator::Greater,
            TokenKind::GreaterEqual => BinaryOperator::GreaterEqual,
            TokenKind::Regex => BinaryOperator::Regex,
            TokenKind::NotRegex => BinaryOperator::NotRegex,
            _ => return Ok(left),
        };
        self.advance()?;
        let right = self.parse_addition()?;
        if matches!(
            self.current.kind,
            TokenKind::Equal
                | TokenKind::NotEqual
                | TokenKind::Less
                | TokenKind::LessEqual
                | TokenKind::Greater
                | TokenKind::GreaterEqual
                | TokenKind::Regex
                | TokenKind::NotRegex
        ) {
            return Err(ParseError::new(
                self.current.offset,
                "comparison operators cannot be chained",
            ));
        }
        Ok(binary(operator, left, right))
    }

    fn parse_addition(&mut self) -> Result<Expression, ParseError> {
        let mut expression = self.parse_multiplication()?;
        loop {
            let operator = match self.current.kind {
                TokenKind::Plus => BinaryOperator::Add,
                TokenKind::Minus => BinaryOperator::Subtract,
                _ => return Ok(expression),
            };
            self.advance()?;
            let right = self.parse_multiplication()?;
            expression = binary(operator, expression, right);
        }
    }

    fn parse_multiplication(&mut self) -> Result<Expression, ParseError> {
        let mut expression = self.parse_unary()?;
        loop {
            let operator = match self.current.kind {
                TokenKind::Star => BinaryOperator::Multiply,
                TokenKind::Slash => BinaryOperator::Divide,
                TokenKind::Percent => BinaryOperator::Modulo,
                _ => return Ok(expression),
            };
            self.advance()?;
            let right = self.parse_unary()?;
            expression = binary(operator, expression, right);
        }
    }

    fn parse_unary(&mut self) -> Result<Expression, ParseError> {
        if self.current.kind == TokenKind::Minus {
            self.advance()?;
            return Ok(Expression::Unary {
                operator: UnaryOperator::Negate,
                expression: Box::new(self.parse_unary()?),
            });
        }
        self.parse_primary()
    }

    fn parse_primary(&mut self) -> Result<Expression, ParseError> {
        let token = self.advance()?;
        match token.kind {
            TokenKind::Number(value) => Ok(Expression::Value(Value::Number(value))),
            TokenKind::String(value) if value == "." => Ok(Expression::Value(Value::Missing)),
            TokenKind::String(value) => Ok(Expression::Value(Value::String(value))),
            TokenKind::File(path) => Ok(Expression::Value(Value::File(path))),
            TokenKind::Identifier(name) if name == "." => Err(ParseError::new(
                token.offset,
                "missing values must be quoted as \".\"",
            )),
            TokenKind::Identifier(name) => self.parse_identifier(name),
            TokenKind::LeftParen => {
                let expression = self.parse_site_or()?;
                self.expect(TokenKind::RightParen, "expected )")?;
                Ok(expression)
            }
            _ => Err(ParseError::new(token.offset, "expected a value")),
        }
    }

    fn parse_identifier(&mut self, mut name: String) -> Result<Expression, ParseError> {
        if self.current.kind == TokenKind::LeftParen {
            self.advance()?;
            let mut arguments = Vec::new();
            if self.current.kind != TokenKind::RightParen {
                loop {
                    arguments.push(self.parse_site_or()?);
                    if self.current.kind != TokenKind::Comma {
                        break;
                    }
                    self.advance()?;
                }
            }
            self.expect(TokenKind::RightParen, "expected ) after function arguments")?;
            return Ok(Expression::Function { name, arguments });
        }

        let namespace = match name.to_ascii_uppercase().as_str() {
            "INFO" if self.current.kind == TokenKind::Slash => Some(Namespace::Info),
            "FMT" | "FORMAT" if self.current.kind == TokenKind::Slash => Some(Namespace::Format),
            _ => None,
        };
        if namespace.is_some() {
            self.advance()?;
            let token = self.advance()?;
            let TokenKind::Identifier(field) = token.kind else {
                return Err(ParseError::new(token.offset, "expected a field name"));
            };
            name = field;
        }
        let subscript = if self.current.kind == TokenKind::LeftBracket {
            Some(self.parse_subscript()?)
        } else {
            None
        };
        Ok(Expression::Value(Value::Field(Field {
            namespace,
            name,
            subscript,
        })))
    }

    fn parse_subscript(&mut self) -> Result<Subscript, ParseError> {
        self.expect(TokenKind::LeftBracket, "expected [")?;
        let first = if self.current.kind == TokenKind::Colon {
            Selector::Any
        } else {
            self.parse_selector()?
        };
        let subscript = if self.current.kind == TokenKind::Colon {
            self.advance()?;
            let second = if self.current.kind == TokenKind::RightBracket {
                Selector::Any
            } else {
                self.parse_selector()?
            };
            Subscript::SampleValues {
                samples: first,
                values: second,
            }
        } else {
            Subscript::Values(first)
        };
        self.expect(TokenKind::RightBracket, "expected ]")?;
        Ok(subscript)
    }

    fn parse_selector(&mut self) -> Result<Selector, ParseError> {
        match self.current.kind.clone() {
            TokenKind::Star => {
                self.advance()?;
                Ok(Selector::Any)
            }
            TokenKind::Identifier(name) if name.eq_ignore_ascii_case("GT") => {
                self.advance()?;
                Ok(Selector::Genotype)
            }
            TokenKind::File(path) => {
                self.advance()?;
                Ok(Selector::File(path))
            }
            TokenKind::Number(_) => self.parse_indices(),
            _ => Err(ParseError::new(
                self.current.offset,
                "expected *, GT, a file, or an index",
            )),
        }
    }

    fn parse_indices(&mut self) -> Result<Selector, ParseError> {
        let mut indices = Vec::new();
        loop {
            let start = self.parse_index()?;
            let end = if self.current.kind == TokenKind::Minus {
                self.advance()?;
                if matches!(self.current.kind, TokenKind::Number(_)) {
                    let end = self.parse_index()?;
                    if end < start {
                        return Err(ParseError::new(
                            self.current.offset,
                            "index range ends before it starts",
                        ));
                    }
                    Some(end)
                } else {
                    Some(usize::MAX)
                }
            } else {
                None
            };
            indices.push(IndexSelector { start, end });
            if self.current.kind != TokenKind::Comma {
                break;
            }
            self.advance()?;
        }
        Ok(Selector::Indices(indices))
    }

    fn parse_index(&mut self) -> Result<usize, ParseError> {
        let token = self.advance()?;
        let TokenKind::Number(value) = token.kind else {
            return Err(ParseError::new(token.offset, "expected an index"));
        };
        if value.fract() != 0.0 || value < 0.0 || value > usize::MAX as f64 {
            return Err(ParseError::new(
                token.offset,
                "index must be a nonnegative integer",
            ));
        }
        Ok(value as usize)
    }

    fn expect(&mut self, expected: TokenKind, message: &str) -> Result<(), ParseError> {
        if self.current.kind != expected {
            return Err(ParseError::new(self.current.offset, message));
        }
        self.advance()?;
        Ok(())
    }
}

fn binary(operator: BinaryOperator, left: Expression, right: Expression) -> Expression {
    Expression::Binary {
        operator,
        left: Box::new(left),
        right: Box::new(right),
    }
}

pub(crate) fn parse(source: &str) -> Result<Expression, ParseError> {
    Parser::new(source)?.parse()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn field(namespace: Option<Namespace>, name: &str) -> Expression {
        Expression::Value(Value::Field(Field {
            namespace,
            name: name.to_owned(),
            subscript: None,
        }))
    }

    fn number(value: f64) -> Expression {
        Expression::Value(Value::Number(value))
    }

    #[test]
    fn arithmetic_precedence_and_unary_negation_are_explicit() {
        assert_eq!(
            parse("QUAL + 2 * -3 >= 1e-4").unwrap(),
            binary(
                BinaryOperator::GreaterEqual,
                binary(
                    BinaryOperator::Add,
                    field(None, "QUAL"),
                    binary(
                        BinaryOperator::Multiply,
                        number(2.0),
                        Expression::Unary {
                            operator: UnaryOperator::Negate,
                            expression: Box::new(number(3.0)),
                        },
                    ),
                ),
                number(0.0001),
            )
        );
    }

    #[test]
    fn sample_and_site_logical_operators_remain_distinct() {
        let expression = parse("FMT/DP>10 & FMT/GQ>20 || QUAL>30 && GT=\"het\"").unwrap();
        let Expression::Binary {
            operator: BinaryOperator::SiteOr,
            left,
            right,
        } = expression
        else {
            panic!("expected site OR");
        };
        assert!(matches!(
            *left,
            Expression::Binary {
                operator: BinaryOperator::SampleAnd,
                ..
            }
        ));
        assert!(matches!(
            *right,
            Expression::Binary {
                operator: BinaryOperator::SiteAnd,
                ..
            }
        ));
    }

    #[test]
    fn namespaces_and_subscripts_cover_samples_values_ranges_and_gt() {
        let expression = parse("FORMAT/AD[0,2-4:GT] > INFO/AF[1-]").unwrap();
        let Expression::Binary { left, right, .. } = expression else {
            panic!("expected comparison");
        };
        let Expression::Value(Value::Field(left)) = *left else {
            panic!("expected FORMAT field");
        };
        assert_eq!(left.namespace, Some(Namespace::Format));
        assert_eq!(left.name, "AD");
        assert_eq!(
            left.subscript,
            Some(Subscript::SampleValues {
                samples: Selector::Indices(vec![
                    IndexSelector {
                        start: 0,
                        end: None,
                    },
                    IndexSelector {
                        start: 2,
                        end: Some(4),
                    },
                ]),
                values: Selector::Genotype,
            })
        );
        let Expression::Value(Value::Field(right)) = *right else {
            panic!("expected INFO field");
        };
        assert_eq!(right.namespace, Some(Namespace::Info));
        assert_eq!(
            right.subscript,
            Some(Subscript::Values(Selector::Indices(vec![IndexSelector {
                start: 1,
                end: Some(usize::MAX),
            }])))
        );
    }

    #[test]
    fn functions_files_wildcards_and_regexes_parse() {
        let expression = parse("phred(binom(FMT/AD[0:*])) > 20 && ID=@ids.txt").unwrap();
        assert!(matches!(
            expression,
            Expression::Binary {
                operator: BinaryOperator::SiteAnd,
                ..
            }
        ));
        let regex = parse("INFO/CSQ ~ \"missense.*deleterious/i\"").unwrap();
        assert!(matches!(
            regex,
            Expression::Binary {
                operator: BinaryOperator::Regex,
                ..
            }
        ));
        assert!(parse("FMT/AD[:1] > 2 & FMT/AD[0:] > 3").is_ok());
    }

    #[test]
    fn missing_values_use_the_bcftools_quoted_spelling() {
        let Expression::Binary { right, .. } = parse("X = \".\"").unwrap() else {
            panic!("expected comparison");
        };
        assert_eq!(*right, Expression::Value(Value::Missing));
        assert!(parse("X = .").is_err());
    }

    #[test]
    fn malformed_expressions_fail_at_the_boundary() {
        for source in [
            "",
            "QUAL >",
            "QUAL > 1 > 0",
            "FMT/AD[] > 1",
            "FMT/AD[3-1] > 1",
            "SUM(DP",
            "ID=@",
            "INFO/ > 1",
            "\"unterminated",
            "1e-",
        ] {
            assert!(parse(source).is_err(), "{source}");
        }
    }
}
