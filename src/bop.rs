use crate::point::Point;

#[allow(non_snake_case)]
pub fn bop(center: &Point, neighbors: &Vec<Point>, bop_params: &[f64], r_cut: f64) -> f64 {
    let A = bop_params[0];
    let a = bop_params[1];
    let alpha = bop_params[2];
    let B = bop_params[3];
    let beta = bop_params[4];
    let h = bop_params[5];
    let lambda = bop_params[6];
    let sigma = bop_params[7];

    let mut E_term_first = 0.0;
    for j in 0..neighbors.len() {
        let r_ij = center.distance(&neighbors[j]);
        let S_ij = calculate_S_ij(center, neighbors, j, r_cut, lambda);
        let z_ij = calcualte_z_ij(center, neighbors, j, r_cut, a, h, lambda);
        let b_ij = (1.0 + z_ij).powf(-0.5);
        let cutoff = cutoff_func(r_ij, r_cut);
        E_term_first += ((A - alpha * r_ij).exp() - S_ij * b_ij * (B - beta * r_ij).exp()) * cutoff;
    }
    E_term_first *= 0.5;
    let E_promo = calculate_E_promo(center, neighbors, r_cut, a, h, sigma, lambda);
    E_term_first + E_promo
}

fn cutoff_func(r: f64, r_cut: f64) -> f64 {
    if r > r_cut {
        0.0
    } else {
        let d_pow = (0.25 * r_cut).powi(4);
        let pow_term = (r - r_cut).powi(4);
        pow_term / (d_pow + pow_term)
    }
}

#[allow(non_snake_case)]
fn calculate_S_ijk(
    center: &Point,
    neighbors: &Vec<Point>,
    j: usize,
    k: usize,
    r_cut: f64,
    lambda: f64,
) -> f64 {
    let r_ij = center.distance(&neighbors[j]);
    let r_ik = center.distance(&neighbors[k]);
    let r_jk = neighbors[j].distance(&neighbors[k]);
    let _temp = r_ik + r_jk - r_ij;
    let S_ijk = 1.0 - cutoff_func(_temp, r_cut) * (-lambda.powi(2) * _temp).exp();
    S_ijk
}

#[allow(non_snake_case)]
fn calculate_S_ij(
    center: &Point,
    neighbors: &Vec<Point>,
    j: usize,
    r_cut: f64,
    lambda: f64,
) -> f64 {
    let mut S_ij = 1.0;
    for k in 0..neighbors.len() {
        if k != j {
            S_ij *= calculate_S_ijk(center, neighbors, j, k, r_cut, lambda);
        }
    }
    S_ij
}

#[allow(non_snake_case)]
fn calcualte_z_ij(
    center: &Point,
    neighbors: &Vec<Point>,
    j: usize,
    r_cut: f64,
    a: f64,
    h: f64,
    lambda: f64,
) -> f64 {
    let mut z_ij = 0.0;
    for k in 0..neighbors.len() {
        if k != j {
            let mut _temp = (center.cos(&neighbors[j], &neighbors[k]) + h).powi(2);
            _temp *= cutoff_func(center.distance(&neighbors[k]), r_cut);
            _temp *= calculate_S_ij(center, neighbors, k, r_cut, lambda);

            z_ij += _temp;
        }
    }
    z_ij *= a.powi(2);
    z_ij
}

#[allow(non_snake_case)]
fn calculate_E_promo(
    center: &Point,
    neighbors: &Vec<Point>,
    r_cut: f64,
    a: f64,
    h: f64,
    sigma: f64,
    lambda: f64,
) -> f64 {
    let mut E_promo = 0.0;
    for j in 0..neighbors.len() {
        let S_ij = calculate_S_ij(center, neighbors, j, r_cut, lambda);
        let z_ij = calcualte_z_ij(center, neighbors, j, r_cut, a, h, lambda);
        let b_ij = (1.0 + z_ij).powf(-0.5);
        let cutoff = cutoff_func(center.distance(&neighbors[j]), r_cut);
        E_promo += b_ij * S_ij * cutoff;
    }
    E_promo = E_promo.sqrt();
    E_promo *= -sigma;
    E_promo
}
