use std::f64::consts::PI;
use num_complex::Complex64;
use ndarray::Array2;
use tightbinding::float::is_near_float;
use tightbinding::tetra::EvecCache;
use tightbinding::dos::{DosValues, dos_from_num};
use model::BlgMoireModel;

pub fn calculate_dos(w: f64, hbar_v: f64, u: f64, t: &Array2<Complex64>, k_max: f64) -> DosValues {
    let model = BlgMoireModel::new(w, hbar_v, u, t.clone());

    let hk_fn = |k| model.hk_lat(&k);

    let dims = [256, 256, 1];

    let k_start = [-k_max, -k_max, 0.0];
    let k_stop = [k_max, k_max, 1.0];

    let num_energies = 1001;

    let cache = EvecCache::new(hk_fn, model.bands(), dims, k_start, k_stop);

    let unscaled_dos = dos_from_num(&cache, num_energies);
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
