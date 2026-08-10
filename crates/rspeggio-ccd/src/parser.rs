// crates/rspeggio-ccd/src/parser.rs
// TODO: remove once load_ccd_file/tokenize_cif_row are called from real (non-test) code,
// and CcdAtom/CcdComponent are actually constructed here.
#![allow(dead_code)]
#![allow(unused_imports)]

use crate::component::{CcdAtom, CcdComponent};

fn load_ccd_file(path: &str) -> Result<String, std::io::Error> {
    std::fs::read_to_string(path)
}

// TODO: in the future this could be zero-copy borrowed return <Vec &str> but for now will copy
// TODO: revisit this as an iterator-based walk (chars.iter().peekable(), .next()/.peek())
// instead of manual index bookkeeping, once quote-handling (step 3) is done and this is
// fully working end to end. Not a correctness issue, just a more idiomatic-Rust refactor
// to come back to deliberately, not mid-implementation.
fn tokenize_cif_row(line: &str) -> Vec<String> {
    let chars: Vec<char> = line.chars().collect();
    let mut tokens: Vec<String> = Vec::new();
    let mut i = 0;
    let mut start: usize = 0;
    while i < chars.len() {
        // 1. skip any whitespace characters, advancing i
        // 2. if we've hit the end of the line, stop
        // 3. otherwise, read characters until the next whitespace (or end of line)
        // 4. collect that range of chars into a String, push onto tokens
        if chars[i] == '"' || chars[i] == '\'' {
            let quote_char = chars[i];
            let mut j = i + 1;
            while j < chars.len() {
                if chars[j] == quote_char && (j + 1 >= chars.len() || chars[j + 1].is_whitespace())
                {
                    break;
                }
                j += 1;
            }
            let token: String = chars[i + 1..j].iter().collect();
            tokens.push(token);
            i = j + 1;
            start = i;
        } else if chars[i].is_whitespace() {
            if start < i {
                let token: String = chars[start..i].iter().collect(); // not inclusive of i
                tokens.push(token);
            }
            start = i + 1;
            i += 1;
        } else {
            i += 1;
        }
    }
    if start < chars.len() {
        let final_token: String = chars[start..chars.len()].iter().collect();
        tokens.push(final_token);
    }
    tokens
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_a_real_ccd_file() {
        let result = load_ccd_file("tests/fixtures/ADP_ideal.cif");
        assert!(result.is_ok());

        let contents = result.expect("fixture file should exist and be readable");
        assert!(contents.contains("ADP"));
    }

    #[test]
    fn empty_line_yields_no_tokens() {
        assert_eq!(tokenize_cif_row(""), Vec::<String>::new());
    }

    #[test]
    fn splits_unquoted_fields_on_whitespace() {
        let tokens = tokenize_cif_row("ALA N N N 0 1 N N N N N N");
        assert_eq!(
            tokens,
            vec!["ALA", "N", "N", "N", "0", "1", "N", "N", "N", "N", "N", "N"]
        );
    }

    #[test]
    fn handles_quoted_fields_with_apostrophes() {
        let tokens = tokenize_cif_row("TPP \"N1'\"  \"N1'\"  N 0 1 Y N N");
        assert_eq!(
            tokens,
            vec!["TPP", "N1'", "N1'", "N", "0", "1", "Y", "N", "N"]
        );
    }
}
