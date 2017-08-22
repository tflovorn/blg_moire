extern crate num_complex;
extern crate ndarray;
extern crate linxal;
#[macro_use]
extern crate serde_json;
extern crate tightbinding;
extern crate blg_moire;

use std::fs::File;
use std::io::Write;
use num_complex::Complex64;
use ndarray::Array2;
use linxal::eigenvalues::SymEigen;
use linxal::types::Symmetric;
use tightbinding::vec_util::transpose_vecs;
use tightbinding::float::is_near_float;
use blg_moire::model::BlgMoireModel;

fn main() {
    // Setting w = 1, we give all energies in units of w.
    let w = 1.0;
    // Additionally setting hbar * v = 1, we give lengths in units of (hbar * v / w).
    let hbar_v = 1.0;

    let us = [0.4, 0.8, 1.2, 10.0];
    let k_maxs = [10.0, 10.0, 10.0, 40.0];
    let ts = BlgMoireModel::ts();

    let num_ks = 1000;

    for (t_index, t) in ts.iter().enumerate() {
        for (u_index, (u, k_max)) in us.iter().zip(k_maxs.iter()).enumerate() {
            let (kxs, kys) = get_ks(num_ks, *k_max);

            let out_prefix_kx = format!("blg_moire_kx_t{}_u{}", t_index, u_index);
            let out_prefix_ky = format!("blg_moire_ky_t{}_u{}", t_index, u_index);

            write_spectrum(w, hbar_v, *u, t_index, t, &kxs, &out_prefix_kx);
            write_spectrum(w, hbar_v, *u, t_index, t, &kys, &out_prefix_ky);
        }
    }
}

fn get_ks(num_ks: usize, k_max: f64) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let k_start = -k_max;
    let k_stop = k_max;
    let step = (k_stop - k_start) / (num_ks as f64);

    let mut kxs = Vec::new();
    let mut kys = Vec::new();

    for i in 0..num_ks + 1 {
        kxs.push([k_start + (i as f64) * step, 0.0, 0.0]);
        kys.push([0.0, k_start + (i as f64) * step, 0.0]);
    }

    (kxs, kys)
}

fn write_spectrum(
    w: f64,
    hbar_v: f64,
    u: f64,
    t_index: usize,
    t: &Array2<Complex64>,
    ks: &Vec<[f64; 3]>,
    out_prefix: &str,
) {
    let model = BlgMoireModel::new(w, hbar_v, u, t.clone());

    let mut ekms = Vec::new();

    for k in ks.iter() {
        let hk = model.hk_lat(k);
        let solution = SymEigen::compute(&hk, Symmetric::Upper, true).unwrap();
        ekms.push(solution.values.to_vec());
    }

    let emks = transpose_vecs(&ekms);

    // Captions assume that w = 1 and hbar * v = 1.
    // Could rescale ks, emks before outputting to avoid this assert.
    assert!(is_near_float(w, 1.0, 1e-12, 1e-12));
    assert!(is_near_float(hbar_v, 1.0, 1e-12, 1e-12));

    let xlabel = r"$k_x \hbar v / w$";
    //let xlabel = r"$k_y \hbar v / w$";
    let ylabel = r"$E / w$";

    // TODO caption for T.
    let caption = format!(r"$T_{} \, ; \, U = {}$", t_index + 1, u);

    let json_out = json!({
        "ks": ks,
        "emks": emks,
        "xlabel": xlabel,
        "ylabel": ylabel,
        "caption": caption,
    });

    let out_path = format!("{}_spectrum.json", out_prefix);

    let mut file = File::create(out_path).expect("Eror creating output file");
    file.write_all(format!("{}", json_out).as_bytes()).expect(
        "Error writing output file",
    );
}
