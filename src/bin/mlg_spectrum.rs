extern crate num_complex;
extern crate ndarray;
extern crate linxal;
#[macro_use]
extern crate serde_json;
extern crate tightbinding;
extern crate blg_moire;

use std::fs::File;
use std::io::Write;
use linxal::eigenvalues::SymEigen;
use linxal::types::Symmetric;
use tightbinding::vec_util::transpose_vecs;
use blg_moire::model::MlgNNModel;

fn main() {
    let t = 1.0;

    let model = MlgNNModel::new(t);

    let num_ks: usize = 1000;
    let k_start = 0.0;
    let k_stop = 1.0;
    let step = (k_stop - k_start) / (num_ks as f64);
    let mut ks = Vec::new();

    for i in 0..num_ks + 1 {
        let k = k_start + (i as f64) * step;
        ks.push([k, k, 0.0]);
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

    let out_path = "mlg_spectrum.json";

    let mut file = File::create(out_path).expect("Eror creating output file");
    file.write_all(format!("{}", json_out).as_bytes()).expect(
        "Error writing output file",
    );
}
