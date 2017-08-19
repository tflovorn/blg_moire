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
use ndarray::arr2;
use linxal::eigenvalues::SymEigen;
use linxal::types::Symmetric;
use tightbinding::vec_util::transpose_vecs;
use blg_moire::model::BlgMoireModel;

fn main() {
    // Setting w = 1, we give all energies in units of w.
    let w = 1.0;
    // Additionally setting hbar * v = 1, we give lengths in units of (hbar * v / w).
    let hbar_v = 1.0;

    // TODO range over U values
    let u = 0.4;

    // TODO range over T values
    let t = arr2(
        &[
            [Complex64::new(-1.0, 0.0), Complex64::new(2.0, 0.0)],
            [Complex64::new(2.0, 0.0), Complex64::new(-1.0, 0.0)],
        ],
    );

    let model = BlgMoireModel::new(w, hbar_v, u, t);

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

    let out_path = "blg_moire_spectrum.json";

    let mut file = File::create(out_path).expect("Eror creating output file");
    file.write_all(format!("{}", json_out).as_bytes()).expect(
        "Error writing output file",
    );
}
