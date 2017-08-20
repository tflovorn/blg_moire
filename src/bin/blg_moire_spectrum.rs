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
use blg_moire::model::BlgMoireModel;

fn main() {
    // Setting w = 1, we give all energies in units of w.
    let w = 1.0;
    // Additionally setting hbar * v = 1, we give lengths in units of (hbar * v / w).
    let hbar_v = 1.0;

    let us = [0.4, 0.8, 1.2];
    let ts = BlgMoireModel::ts();

    for (t_index, t) in ts.iter().enumerate() {
        for (u_index, u) in us.iter().enumerate() {
            let out_path = format!("blg_moire_t{}_u{}_spectrum.json", t_index, u_index);

            write_spectrum(w, hbar_v, *u, t, &out_path);
        }
    }
}

fn write_spectrum(w: f64, hbar_v: f64, u: f64, t: &Array2<Complex64>, out_path: &str) {
    let model = BlgMoireModel::new(w, hbar_v, u, t.clone());

    let num_ks: usize = 1000;
    let k_start = -10.0;
    let k_stop = 10.0;
    let step = (k_stop - k_start) / (num_ks as f64);
    let mut ks = Vec::new();

    for i in 0..num_ks + 1 {
        ks.push([k_start + (i as f64) * step, 0.0, 0.0]);
    }

    let mut ekms = Vec::new();

    for k in ks.iter() {
        let hk = model.hk_lat(k);
        let solution = SymEigen::compute(&hk, Symmetric::Upper, true).unwrap();
        ekms.push(solution.values.to_vec());
    }

    let emks = transpose_vecs(&ekms);

    let json_out = json!({
        "ks": ks,
        "emks": emks,
    });

    let mut file = File::create(out_path).expect("Eror creating output file");
    file.write_all(format!("{}", json_out).as_bytes()).expect(
        "Error writing output file",
    );
}
