use std::f64::consts::PI;
use num_complex::Complex64;
use ndarray::{Array2, arr2};

#[derive(Clone)]
pub struct BlgMoireModel {
    w: f64,
    hbar_v: f64,
    u: f64,
    t: Array2<Complex64>,
    sigma_t: Vec<Array2<Complex64>>,
    sigma_b: Vec<Array2<Complex64>>,
}

impl BlgMoireModel {
    pub fn new(theta: f64, w: f64, hbar_v: f64, u: f64, t: Array2<Complex64>) -> BlgMoireModel {
        BlgMoireModel {
            w,
            hbar_v,
            u,
            t,
            sigma_t: pauli_matrices(-theta / 2.0),
            sigma_b: pauli_matrices(theta / 2.0),
        }
    }

    pub fn ts() -> Vec<Array2<Complex64>> {
        let t0 = arr2(
            &[
                [Complex64::new(-1.0, 0.0), Complex64::new(2.0, 0.0)],
                [Complex64::new(2.0, 0.0), Complex64::new(-1.0, 0.0)],
            ],
        );

        let t1 = arr2(
            &[
                [Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)],
                [Complex64::new(3.0, 0.0), Complex64::new(0.0, 0.0)],
            ],
        );

        let rt3 = 3.0_f64.sqrt();
        let t2 = arr2(
            &[
                [Complex64::new(1.0 - rt3, 0.0), Complex64::new(1.0, 0.0)],
                [
                    Complex64::new(1.0 + rt3, 0.0),
                    Complex64::new(1.0 - rt3, 0.0),
                ],
            ],
        );

        vec![t0, t1, t2]
    }

    pub fn t(x: f64, y: f64) -> Array2<Complex64> {
        let e4 = Complex64::new(0.0, 4.0 * PI * y / 3.0).exp();
        let e2 = Complex64::new(0.0, -2.0 * PI * y / 3.0).exp();
        let xval = 2.0 * PI * x / 3.0_f64.sqrt();

        let td = e4 - 2.0 * e2 * xval.cos();
        let tba = e4 + 2.0 * e2 * (xval + PI / 6.0).sin();
        let tab = e4 + 2.0 * e2 * (xval + PI / 3.0).cos();

        arr2(&[[td, tab], [tba, td]])
    }

    pub fn hk_lat(&self, k_lat: &[f64; 3]) -> Array2<Complex64> {
        let k_dot_sigma_t = (&self.sigma_t[0] * k_lat[0] + &self.sigma_t[1] * k_lat[1]) *
            self.hbar_v;
        let k_dot_sigma_b = (&self.sigma_b[0] * k_lat[0] + &self.sigma_b[1] * k_lat[1]) *
            self.hbar_v;
        let ui = Array2::eye(2) * Complex64::new(self.u, 0.0) * self.w;

        let mut hk = Array2::zeros((4, 4));

        hk.slice_mut(s![0..2, 0..2]).assign(&(&k_dot_sigma_t - &ui));
        hk.slice_mut(s![0..2, 2..4]).assign(&(&self.t * self.w));
        hk.slice_mut(s![2..4, 0..2]).assign(
            &((&self.t * self.w).t().map(
                |x| {
                    x.conj()
                },
            )),
        );
        hk.slice_mut(s![2..4, 2..4]).assign(&(&k_dot_sigma_b + &ui));

        hk
    }

    pub fn bands(&self) -> usize {
        4
    }
}

pub fn pauli_matrices(theta: f64) -> Vec<Array2<Complex64>> {
    let sigma_x = arr2(
        &[
            [Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
            [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
        ],
    );
    let sigma_y = arr2(
        &[
            [Complex64::new(0.0, 0.0), Complex64::new(0.0, -1.0)],
            [Complex64::new(0.0, 1.0), Complex64::new(0.0, 0.0)],
        ],
    );
    let sigma_z = arr2(
        &[
            [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
            [Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0)],
        ],
    );

    let sigma_x_rot = &sigma_x * (theta / 2.0).cos() + &sigma_y * (theta / 2.0).sin();
    let sigma_y_rot = &sigma_y * (theta / 2.0).cos() - &sigma_x * (theta / 2.0).sin();

    vec![sigma_x_rot, sigma_y_rot, sigma_z]
}

#[cfg(test)]
mod tests {
    use super::BlgMoireModel;
    use tightbinding::float::is_near_complex;

    #[test]
    fn check_t() {
        let ts_fixed = BlgMoireModel::ts();
        let ts_calc = vec![
            BlgMoireModel::t(0.0, 0.0),
            BlgMoireModel::t(1.0 / (2.0 * 3.0_f64.sqrt()), 0.0),
            BlgMoireModel::t(1.0 / (4.0 * 3.0_f64.sqrt()), 0.0),
        ];

        for (t_fixed, t_calc) in ts_fixed.iter().zip(ts_calc.iter()) {
            for (t_fixed_elem, t_calc_elem) in t_fixed.iter().zip(t_calc.iter()) {
                assert!(is_near_complex(*t_fixed_elem, *t_calc_elem, 1e-12, 1e-12));
            }
        }
    }
}
