extern crate num_complex;
extern crate ndarray;
#[macro_use]
extern crate serde_json;
extern crate tightbinding;
extern crate blg_moire;

use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use num_complex::Complex64;
use ndarray::Array2;
use tightbinding::float::is_near_float;
use tightbinding::tetra::EvecCache;
use tightbinding::dos::{DosValues, dos_from_num};
use blg_moire::model::BlgMoireModel;

fn main() {
    // Setting w = 1, we give all energies in units of w.
    let w = 1.0;
    // Additionally setting hbar * v = 1, we give lengths in units of (hbar * v / w).
    let hbar_v = 1.0;

    let us = [0.4, 0.8, 1.2, 10.0];
    let k_maxs = [10.0, 10.0, 10.0, 40.0];
    let ts = BlgMoireModel::ts();

    for (t_index, t) in ts.iter().enumerate() {
        for (u_index, (u, k_max)) in us.iter().zip(k_maxs.iter()).enumerate() {
            let out_path = format!("blg_moire_t{}_u{}_dos.json", t_index, u_index);

            write_dos(w, hbar_v, *u, t_index, t, *k_max, &out_path);
        }
    }
}

fn write_dos(
    w: f64,
    hbar_v: f64,
    u: f64,
    t_index: usize,
    t: &Array2<Complex64>,
    k_max: f64,
    out_path: &str,
) {
    let model = BlgMoireModel::new(w, hbar_v, u, t.clone());

    let hk_fn = |k| model.hk_lat(&k);

    let dims = [256, 256, 1];

    let k_start = [-k_max, -k_max, 0.0];
    let k_stop = [k_max, k_max, 1.0];

    // Scaling of DOS assumes k_stop[2] - k_start[2] = 1.0.
    assert!(is_near_float(k_start[2], 0.0, 1e-12, 1e-12));
    assert!(is_near_float(k_stop[2], 1.0, 1e-12, 1e-12));

    let num_energies = 1001;

    let cache = EvecCache::new(hk_fn, model.bands(), dims, k_start, k_stop);

    let dos = dos_from_num(&cache, num_energies);
    let dos = rescale_dos(dos, &k_start, &k_stop);

    // Scaling given in captions assumes that w = 1 and hbar * v = 1.
    // Could rescale ks, emks, DOS before outputting to avoid this assert.
    assert!(is_near_float(w, 1.0, 1e-12, 1e-12));
    assert!(is_near_float(hbar_v, 1.0, 1e-12, 1e-12));

    let xlabel = r"$E / w$";
    let ylabel = r"DOS [$(\frac{w}{\hbar v})^2 / w$]";

    let caption = format!(r"$T_{} \, ; \, U / w = {} $", t_index + 1, u);

    let json_out = json!({
        "es": dos.es,
        "total_dos": dos.total_dos,
        "orbital_dos": dos.orbital_dos,
        "xlabel": xlabel,
        "ylabel": ylabel,
        "caption": caption,
    });

    let mut file = File::create(out_path).expect("Eror creating output file");
    file.write_all(format!("{}", json_out).as_bytes()).expect(
        "Error writing output file",
    );
}

fn rescale_dos(dos: DosValues, k_start: &[f64; 3], k_stop: &[f64; 3]) -> DosValues {
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
