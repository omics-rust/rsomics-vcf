use super::{Field, Number, Type};

pub(crate) fn info_definition(key: &str, version: Option<(u8, u8)>) -> Option<Field> {
    let field = match key {
        "AA" => Field::new(Number::Count(1), Type::String),
        "AC" => Field::new(Number::Alternate, Type::Integer),
        "AD" | "ADF" | "ADR" => Field::new(Number::Alleles, Type::Integer),
        "AF" => Field::new(Number::Alternate, Type::Float),
        "AN" | "DP" | "MQ0" | "NS" | "END" => Field::new(Number::Count(1), Type::Integer),
        "BQ" | "MQ" => Field::new(Number::Count(1), Type::Float),
        "CIGAR" => Field::new(Number::Alternate, Type::String),
        "DB" | "H2" | "H3" | "SOMATIC" | "VALIDATED" | "1000G" | "IMPRECISE" | "NOVEL" => {
            Field::new(Number::Count(0), Type::Flag)
        }
        "SB" => Field::new(Number::Count(1), Type::Float),
        "SVTYPE" => Field::new(Number::Count(1), Type::String),
        "SVLEN" if version.is_some_and(|version| version >= (4, 4)) => {
            Field::new(Number::Alternate, Type::Integer)
        }
        "SVLEN" => Field::new(Number::Count(1), Type::Integer),
        "CIPOS" | "CIEND" => Field::new(Number::Count(2), Type::Integer),
        _ => return None,
    };
    Some(field)
}

pub(crate) fn format_definition(key: &str, version: Option<(u8, u8)>) -> Option<Field> {
    let field = match key {
        "AD" | "ADF" | "ADR" => Field::new(Number::Alleles, Type::Integer),
        "DP" | "GQ" | "MQ" | "PQ" | "PS" => Field::new(Number::Count(1), Type::Integer),
        "EC" => Field::new(Number::Alternate, Type::Integer),
        "FT" => Field::new(Number::Count(1), Type::String),
        "GL" | "GP" => Field::new(Number::Genotypes, Type::Float),
        "GLE" => Field::new(Number::Genotypes, Type::String),
        "GT" => Field::new(Number::Count(1), Type::String),
        "HQ" => Field::new(Number::Count(2), Type::Integer),
        "PL" | "PP" => Field::new(Number::Genotypes, Type::Integer),
        "PSL" if version.is_some_and(|version| version >= (4, 4)) => {
            Field::new(Number::Ploidy, Type::String)
        }
        "PSO" | "PSQ" if version.is_some_and(|version| version >= (4, 4)) => {
            Field::new(Number::Ploidy, Type::Integer)
        }
        _ => return None,
    };
    Some(field)
}
