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

    pub fn hk_lat(&self, k_lat: &[f64; 3]) -> Array2<Complex64> {
        let k_dot_sigma = (&self.sigma[0] * k_lat[0] + &self.sigma[1] * k_lat[1]) * self.hbar_v;
        let ui = Array2::eye(2) * Complex64::new(self.u, 0.0) * self.w;

        let mut hk = Array2::zeros((4, 4));

        hk.slice_mut(s![0..2, 0..2]).assign(&(&k_dot_sigma - &ui));
        hk.slice_mut(s![0..2, 2..4]).assign(&(&self.t * self.w));
        hk.slice_mut(s![2..4, 0..2]).assign(
            &((&self.t * self.w).t().map(
                |x| {
                    x.conj()
                },
            )),
        );
        hk.slice_mut(s![2..4, 2..4]).assign(&(&k_dot_sigma + &ui));

        hk
    }

    pub fn bands(&self) -> usize {
        4
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
