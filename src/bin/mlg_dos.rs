extern crate num_complex;
extern crate ndarray;
#[macro_use]
extern crate serde_json;
extern crate tightbinding;
extern crate blg_moire;

use std::fs::File;
use std::io::Write;
use tightbinding::tetra::EvecCache;
use tightbinding::dos::dos_from_num;
use blg_moire::model::MlgNNModel;

fn main() {
    let t = 1.0;

    let bands = 2;

    let model = MlgNNModel::new(t);

    let hk_fn = |k| model.hk_lat(&k);

    let use_curvature_correction = true;
    let dims = [256, 256, 1];
    let k_start = [0.0, 0.0, 0.0];
    let k_stop = [1.0, 1.0, 1.0];
    let num_energies = 1000;

    let cache = EvecCache::new(hk_fn, bands, dims, k_start, k_stop);

    let (es, dos) = dos_from_num(&cache, num_energies, use_curvature_correction);

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

    let out_path = "mlg_dos.json";

    let mut file = File::create(out_path).expect("Eror creating output file");
    file.write_all(format!("{}", json_out).as_bytes()).expect(
        "Error writing output file",
    );
}
