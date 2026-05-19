use log::info;
use pyo3::prelude::*;

#[pyfunction]
pub fn sum_as_string(a: &str, b: &str) -> String {
    info!("Adding the numbers - {} and {}", a, b);
    a.to_owned() + b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sum_as_string() {
        assert_eq!(sum_as_string("1", "2"), "12");
        assert_eq!(sum_as_string("foo", "bar"), "foobar");
        assert_eq!(sum_as_string("", "x"), "x");
    }
}
