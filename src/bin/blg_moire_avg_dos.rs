extern crate num_complex;
extern crate ndarray;
#[macro_use]
extern crate serde_json;
extern crate tightbinding;
extern crate blg_moire;

use std::fs::File;
use std::io::Write;
use tightbinding::float::is_near_float;
use tightbinding::dos::DosValues;
use blg_moire::model::BlgMoireModel;
use blg_moire::dos::{calculate_dos, average_dos};

fn main() {
    // Setting w = 1, we give all energies in units of w.
    let w = 1.0;
    // Additionally setting hbar * v = 1, we give lengths in units of (hbar * v / w).
    let hbar_v = 1.0;

    // Scaling given in captions assumes that w = 1 and hbar * v = 1.
    // Could rescale ks, emks, DOS before outputting to avoid this assert.
    assert!(is_near_float(w, 1.0, 1e-12, 1e-12));
    assert!(is_near_float(hbar_v, 1.0, 1e-12, 1e-12));

    let us = [0.1];
    let k_maxs = [10.0];

    // TODO - possible to re-use work from smaller r0's in larger r0's if
    // xy values sampled at larger r0 are strict superset of those at smaller r0.
    // Simple way to implement this: specify n_sample only for largest r0,
    // filter results to get smaller r0's.
    let n_samples = [4, 8, 12];
    let r0s = [0.05, 0.10, 0.15];

    for (u_index, (u, k_max)) in us.iter().zip(k_maxs.iter()).enumerate() {
        for (r0_index, (r0, n_sample)) in r0s.iter().zip(n_samples.iter()).enumerate() {
            let xys = sample_xys(*n_sample, *r0);

            let mut all_dos = Vec::with_capacity(xys.len());

            for (xy_index, &(x, y)) in xys.iter().enumerate() {
                println!(
                    "Calculating DOS for (x, y) = ({}, {}); sample {} of {}",
                    x,
                    y,
                    xy_index + 1,
                    xys.len()
                );
                let t = BlgMoireModel::t(x, y);
                let dos = calculate_dos(w, hbar_v, *u, &t, *k_max);

                let caption = format!(r"$(x, y) = ({:.4}, {:.4}) \, ; \, U / w = {} $", x, y, u);
                let out_path = format!(
                    "blg_moire_val_r0_{}_xy_{}_u_{}_dos.json",
                    r0_index,
                    xy_index,
                    u_index
                );

                write_dos(&dos, &caption, &out_path);

                all_dos.push(dos);
            }

            let avg_dos = average_dos(&all_dos);

            let caption = format!(r"$r_0 = {} \, ; \, U / w = {} $", r0, u);
            let out_path = format!("blg_moire_avg_r0_{}_u_{}_dos.json", r0_index, u_index);

            write_dos(&avg_dos, &caption, &out_path);
        }
    }
}

fn sample_xys(n: usize, r0: f64) -> Vec<(f64, f64)> {
    let step = r0 / (n as f64);

    let mut all_xys = Vec::with_capacity((2 * n + 1).pow(2));

    for j in 0..2 * n + 1 {
        for i in 0..2 * n + 1 {
            all_xys.push((-r0 + (i as f64) * step, -r0 + (j as f64) * step));
        }
    }

    let eps_abs = 1e-12;
    all_xys
        .into_iter()
        .filter(|&(x, y)| (x.powi(2) + y.powi(2)).sqrt() - r0 <= eps_abs)
        .collect()
}

fn write_dos(dos: &DosValues, caption: &str, out_path: &str) {
    let xlabel = r"$E / w$";
    let ylabel = r"DOS [$(\frac{w}{\hbar v})^2 / w$]";

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
