use std::f64::consts::PI;
use num_complex::Complex64;
use ndarray::{Array2, arr2};

#[derive(Clone)]
pub struct BlgMoireModel {
    w: f64,
    hbar_v: f64,
    u: f64,
    t: Array2<Complex64>,
    sigma: Vec<Array2<Complex64>>,
}

impl BlgMoireModel {
    pub fn new(w: f64, hbar_v: f64, u: f64, t: Array2<Complex64>) -> BlgMoireModel {
        BlgMoireModel {
            w,
            hbar_v,
            u,
            t,
            sigma: pauli_matrices(),
        }
    }

    pub fn hk_lat(&self, k_lat: &[f64; 3]) -> Array2<Complex64> {
        let k_dot_sigma = &self.sigma[0] * k_lat[0] + &self.sigma[1] * k_lat[1];
        let ui = Array2::eye(2) * Complex64::new(self.u, 0.0);

        let mut hk = Array2::zeros((4, 4));

        hk.slice_mut(s![0..2, 0..2]).assign(&(&k_dot_sigma - &ui));
        hk.slice_mut(s![0..2, 2..4]).assign(&self.t);
        hk.slice_mut(s![2..4, 0..2]).assign(
            &(self.t.t().map(|x| x.conj())),
        );
        hk.slice_mut(s![2..4, 2..4]).assign(&(&k_dot_sigma + &ui));

        hk
    }
}

pub fn pauli_matrices() -> Vec<Array2<Complex64>> {
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

    vec![sigma_x, sigma_y, sigma_z]
}

#[derive(Clone)]
pub struct MlgNNModel {
    t: f64,
}

impl MlgNNModel {
    pub fn new(t: f64) -> MlgNNModel {
        MlgNNModel { t }
    }

    pub fn hk_lat(&self, k_lat: &[f64; 3]) -> Array2<Complex64> {
        let mut hk = Array2::zeros((2, 2));

        hk[[0, 1]] = 1.0 + (Complex64::i() * 2.0 * PI * k_lat[1]).exp() +
            (Complex64::i() * 2.0 * PI * (k_lat[0] + k_lat[1])).exp();
        hk[[1, 0]] = hk[[0, 1]].conj();

        hk
    }
}
