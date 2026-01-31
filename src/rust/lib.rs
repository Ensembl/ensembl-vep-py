use pyo3::prelude::*;

// plain Rust helper, not exposed to Python: useful for threading and unit test
fn sum_as_string_impl(a: &str, b: &str) -> String {
    a.to_owned() + b
}

#[pymodule]
mod _ensembl_vep_core {
    use super::sum_as_string_impl;
    use pyo3::prelude::*;

    /// Can we do this more elegantly such as PyModule.add_function?
    #[pyfunction]
    fn sum_as_string(a: &str, b: &str) -> PyResult<String> {
        Ok(sum_as_string_impl(a, b))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sum_as_string_impl() {
        assert_eq!(sum_as_string_impl("1", "2"), "12");
        assert_eq!(sum_as_string_impl("foo", "bar"), "foobar");
        assert_eq!(sum_as_string_impl("", "x"), "x");
    }
}
