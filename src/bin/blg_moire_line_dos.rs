#[macro_use]
extern crate serde_json;
extern crate tightbinding;
extern crate blg_moire;

use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use tightbinding::float::is_near_float;
use tightbinding::dos::DosValues;
use blg_moire::model::BlgMoireModel;
use blg_moire::dos::calculate_dos;

fn main() {
    // Setting w = 1, we give all energies in units of w.
    let w = 1.0;
    // Additionally setting hbar * v = 1, we give lengths in units of (hbar * v / w).
    let hbar_v = 1.0;
    // Rotation angle between layers.
    let theta = 0.245 * PI / 180.0;

    // Scaling given in captions assumes that w = 1 and hbar * v = 1.
    // Could rescale ks, emks, DOS before outputting to avoid this assert.
    assert!(is_near_float(w, 1.0, 1e-12, 1e-12));
    assert!(is_near_float(hbar_v, 1.0, 1e-12, 1e-12));

    let us = [0.24];
    let k_maxs = [10.0];

    for (u_index, (u, k_max)) in us.iter().zip(k_maxs.iter()).enumerate() {
        let xys = get_xy_list();

        for (xy_index, &(x, y)) in xys.iter().enumerate() {
            println!(
                "Calculating DOS for (x, y) = ({}, {}); sample {} of {}",
                x,
                y,
                xy_index + 1,
                xys.len()
            );
            let t = BlgMoireModel::t(x, y);
            let energy_bounds = Some((-5.0, 5.0));
            let dos = calculate_dos(theta, w, hbar_v, *u, &t, *k_max, energy_bounds);

            let caption = format!(r"$(x, y) = ({:.4}, {:.4}) \, ; \, U / w = {} $", x, y, u);
            let out_path = format!("blg_moire_val_xy_{}_u_{}_dos.json", xy_index, u_index);

            write_dos(&dos, &caption, &out_path);
        }
    }
}

fn get_xy_list() -> Vec<(f64, f64)> {
    let nys = 11;
    let mut xys = Vec::with_capacity(nys);

    let x = 0.0;
    let y_max = 0.25;
    let y_step = y_max / ((nys as f64) - 1.0);

    for iy in 0..nys {
        let y = (iy as f64) * y_step;
        xys.push((x, y));
    }

    xys
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
