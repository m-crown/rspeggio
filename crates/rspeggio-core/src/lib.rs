pub mod config;
pub mod contacts;
pub mod export;
pub mod features;
pub mod hydrogenate;
pub mod join;
pub mod perception;
pub mod rings;
pub mod selection;
pub mod typing;
pub mod unsatisfied;

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
