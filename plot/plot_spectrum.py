import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt

def find_with_prefix(prefix):
    files = [str(x.name) for x in Path('.').iterdir() if x.is_file()]

    with_prefix = []
    for f in files:
        if f.startswith(prefix) and f.endswith("_spectrum.json"):
            full_prefix = f.rpartition("_")[0]
            with_prefix.append((f, full_prefix))

    return with_prefix

def _main():
    parser = argparse.ArgumentParser("plot specta")
    parser.add_argument("prefix", type=str, help="Prefix for calculation")
    parser.add_argument("k_axis", type=int, help="Axis of k-point distribution (0 or 1)")
    args = parser.parse_args()

    paths = find_with_prefix(args.prefix)

    for bands_path, full_prefix in paths:
        with open(bands_path, 'r') as fp:
            bands_data = json.load(fp)

        ks = [k[args.k_axis] for k in bands_data["ks"]]

        if "xlabel" in bands_data:
            plt.xlabel(bands_data["xlabel"])

        if "ylabel" in bands_data:
            plt.ylabel(bands_data["ylabel"])

        if "caption" in bands_data:
            plt.title(bands_data["caption"])

        for em in bands_data["emks"]:
            plt.plot(ks, em, 'k-')

        plot_path = "{}_spectrum.png".format(full_prefix)

        plt.savefig(plot_path, dpi=500, bbox_inches='tight')
        plt.clf()

if __name__ == "__main__":
    _main()
