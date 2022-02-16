pub mod bop;
pub mod point;
use crate::bop::bop;
use crate::point::Point;
use pyo3::prelude::*;

/// Wrapper for crate::bop::bop.
#[pyfunction]
fn calculate_bop(
    center: Vec<f64>,
    neighbors: Vec<Vec<f64>>,
    bop_params: Vec<f64>,
    r_cut: f64,
) -> f64 {
    let center = Point::from_vec(center);
    let neighbors = neighbors
        .into_iter()
        .map(|p| Point::from_vec(p))
        .collect::<Vec<Point>>();
    let energy = bop(&center, &neighbors, &bop_params, r_cut);
    energy
}

#[pymodule]
fn pybop(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calculate_bop, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
