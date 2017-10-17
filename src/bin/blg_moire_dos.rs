extern crate num_complex;
extern crate ndarray;
extern crate itertools;
#[macro_use]
extern crate serde_json;
extern crate tightbinding;
extern crate blg_moire;

use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use num_complex::Complex64;
use ndarray::{Array, Array2};
use itertools::multizip;
use tightbinding::float::is_near_float;
use blg_moire::model::BlgMoireModel;
use blg_moire::dos::calculate_dos;

fn main() {
    // Setting w = 1, we give all energies in units of w.
    let w = 1.0;
    // Additionally setting hbar * v = 1, we give lengths in units of (hbar * v / w).
    let hbar_v = 1.0;
    // Rotation angle between layers.
    //let theta = 0.245 * PI / 180.0;
    let theta = 0.0;

    //let us = [0.1, 0.4, 0.8, 1.2, 10.0];
    //let k_maxs = [10.0, 10.0, 10.0, 10.0, 40.0];
    let us = Array::linspace(0.0, 1.0, 11);
    let k_maxs = [10.0; 11];
    let all_energy_bounds = [Some((-5.0, 5.0)); 11];

    let all_ts = BlgMoireModel::ts();
    let ts = vec![all_ts[0].clone(), all_ts[1].clone()];

    for (t_index, t) in ts.iter().enumerate() {
        for (u_index, (u, k_max, energy_bounds)) in
            multizip((us.iter(), k_maxs.iter(), all_energy_bounds.iter())).enumerate()
        {
            let out_path = format!("blg_moire_t{}_u{}_dos.json", t_index, u_index);

            write_dos(
                theta,
                w,
                hbar_v,
                *u,
                t_index,
                t,
                *k_max,
                *energy_bounds,
                &out_path,
            );
        }
    }
}

fn write_dos(
    theta: f64,
    w: f64,
    hbar_v: f64,
    u: f64,
    t_index: usize,
    t: &Array2<Complex64>,
    k_max: f64,
    energy_bounds: Option<(f64, f64)>,
    out_path: &str,
) {
    let dos = calculate_dos(theta, w, hbar_v, u, t, k_max, energy_bounds);

    // Scaling given in captions assumes that w = 1 and hbar * v = 1.
    // Could rescale ks, emks, DOS before outputting to avoid this assert.
    assert!(is_near_float(w, 1.0, 1e-12, 1e-12));
    assert!(is_near_float(hbar_v, 1.0, 1e-12, 1e-12));

    let xlabel = r"$E / w$";
    let ylabel = r"DOS [$(\frac{w}{\hbar v})^2 / w$]";

    let theta_deg = theta * 180.0 / PI;
    let caption = format!(
        r"$T_{} \, ; \, U / w = {:.2} \, ; \, \theta = {:.3}$ deg.",
        t_index + 1,
        u,
        theta_deg
    );

    let json_out = json!({
        "es": dos.es,
        "total_dos": dos.total_dos,
        "orbital_dos": dos.orbital_dos,
        "xlabel": xlabel,
        "ylabel": ylabel,
        "caption": caption,
        "t_index": t_index,
        "u": u,
        "theta_deg": theta_deg,
    });

    let mut file = File::create(out_path).expect("Eror creating output file");
    file.write_all(format!("{}", json_out).as_bytes()).expect(
        "Error writing output file",
    );
}
