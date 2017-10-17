import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt

def find_with_prefix(prefix):
    files = [str(x.name) for x in Path('.').iterdir() if x.is_file()]

    with_prefix = []
    for f in files:
        if f.startswith(prefix) and f.endswith("_dos.json"):
            full_prefix = f.rpartition("_")[0]
            with_prefix.append((f, full_prefix))

    return with_prefix

def _main():
    parser = argparse.ArgumentParser("Plot total DOS")
    parser.add_argument("prefix", type=str, help="Calculation prefix")
    args = parser.parse_args()

    paths = find_with_prefix(args.prefix)

    for dos_path, full_prefix in paths:
        with open(dos_path) as fp:
            dos_data = json.load(fp)

        es = dos_data["es"]
        total_dos = dos_data["total_dos"]

        if "xlabel" in dos_data:
            plt.xlabel(dos_data["xlabel"])

        if "ylabel" in dos_data:
            plt.ylabel(dos_data["ylabel"])

        if "caption" in dos_data:
            plt.title(dos_data["caption"])

        plot_path = "{}_dos.png".format(full_prefix)

        plt.xlim(min(es), max(es))
        plt.ylim(0.0, max(total_dos))

        plt.plot(es, total_dos, 'k-')

        plt.savefig(plot_path, dpi=500, bbox_inches='tight')
        plt.clf()

if __name__ == "__main__":
    _main()
