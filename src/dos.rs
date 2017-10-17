use std::f64::consts::PI;
use num_complex::Complex64;
use ndarray::Array2;
use tightbinding::float::is_near_float;
use tightbinding::tetra::EvecCache;
use tightbinding::dos::{DosValues, dos_from_num};
use model::BlgMoireModel;

pub fn calculate_dos(
    theta: f64,
    w: f64,
    hbar_v: f64,
    u: f64,
    t: &Array2<Complex64>,
    k_max: f64,
    energy_bounds: Option<(f64, f64)>,
) -> DosValues {
    let model = BlgMoireModel::new(theta, w, hbar_v, u, t.clone());

    let hk_fn = |k| model.hk_lat(&k);

    let dims = [256, 256, 1];

    let k_start = [-k_max, -k_max, 0.0];
    let k_stop = [k_max, k_max, 1.0];

    let num_energies = 1001;

    let cache = EvecCache::new(hk_fn, model.bands(), dims, k_start, k_stop);

    let unscaled_dos = dos_from_num(&cache, num_energies, energy_bounds);
    rescale_dos(unscaled_dos, &k_start, &k_stop)
}

fn rescale_dos(dos: DosValues, k_start: &[f64; 3], k_stop: &[f64; 3]) -> DosValues {
    // Scaling of DOS assumes k_stop[2] - k_start[2] = 1.0.
    assert!(is_near_float(k_start[2], 0.0, 1e-12, 1e-12));
    assert!(is_near_float(k_stop[2], 1.0, 1e-12, 1e-12));

    let scale_factor = (k_stop[0] - k_start[0]) * (k_stop[1] - k_start[1]) / (2.0 * PI).powi(2);
    let rescaled_total = dos.total_dos.iter().map(|x| x * scale_factor).collect();
    let rescaled_orbital = dos.orbital_dos
        .iter()
        .map(|orb| orb.iter().map(|x| x * scale_factor).collect())
        .collect();

    DosValues {
        es: dos.es.clone(),
        total_dos: rescaled_total,
        orbital_dos: rescaled_orbital,
    }
}

pub fn average_dos(all_dos: &Vec<DosValues>) -> DosValues {
    assert!(all_dos.len() > 0);

    let num_es = all_dos[0].es.len();
    let num_orbitals = all_dos[0].orbital_dos.len();
    let num_dos = all_dos.len() as f64;

    let mut total_dos = vec![0.0; num_es];
    let mut orbital_dos = Vec::with_capacity(num_orbitals);

    for _ in 0..num_orbitals {
        orbital_dos.push(vec![0.0; num_es]);
    }

    for dos in all_dos {
        for (e_index, dos_value) in dos.total_dos.iter().enumerate() {
            total_dos[e_index] += dos_value / num_dos;
        }

        for (orbital_index, orbital_dos_values) in dos.orbital_dos.iter().enumerate() {
            for (e_index, dos_value) in orbital_dos_values.iter().enumerate() {
                orbital_dos[orbital_index][e_index] += dos_value / num_dos;
            }
        }
    }

    DosValues {
        es: all_dos[0].es.clone(),
        total_dos,
        orbital_dos,
    }
}
