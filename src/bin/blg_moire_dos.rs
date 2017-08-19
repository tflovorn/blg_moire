extern crate num_complex;
extern crate ndarray;
#[macro_use]
extern crate serde_json;
extern crate tightbinding;
extern crate blg_moire;

use std::fs::File;
use std::io::Write;
use num_complex::Complex64;
use ndarray::arr2;
use tightbinding::tetra::EvecCache;
use tightbinding::dos::dos_from_num;
use blg_moire::model::BlgMoireModel;

fn main() {
    // Setting w = 1, we give all energies in units of w.
    let w = 1.0;
    // Additionally setting hbar * v = 1, we give lengths in units of (hbar * v / w).
    let hbar_v = 1.0;

    // TODO range over U values
    let u = 0.4;

    let bands = 4;

    // TODO range over T values
    let t = arr2(
        &[
            [Complex64::new(-1.0, 0.0), Complex64::new(2.0, 0.0)],
            [Complex64::new(2.0, 0.0), Complex64::new(-1.0, 0.0)],
        ],
    );

    let model = BlgMoireModel::new(w, hbar_v, u, t);

    let hk_fn = |k| model.hk_lat(&k);

    let dims = [128, 128, 1];
    let k_start = [-10.0, -10.0, -10.0];
    let k_stop = [10.0, 10.0, 10.0];
    let num_energies = 1000;

    let cache = EvecCache::new(hk_fn, bands, dims, k_start, k_stop);

    let (es, dos) = dos_from_num(&cache, num_energies);

    let mut total_dos = vec![0.0; es.len()];
    for band_dos in dos.iter() {
        for (e_index, e_dos) in band_dos.iter().enumerate() {
            total_dos[e_index] += *e_dos;
        }
    }

    let json_out = json!({
        "es": es,
        "total_dos": total_dos,
    });

    let out_path = "blg_moire_dos.json";

    let mut file = File::create(out_path).expect("Eror creating output file");
    file.write_all(format!("{}", json_out).as_bytes()).expect(
        "Error writing output file",
    );
}
