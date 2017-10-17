import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt

def find_with_prefix(prefix, dir_path):
    files = sorted([x for x in Path(dir_path).iterdir() if x.is_file()])

    with_prefix = []
    for f in files:
        if f.name.startswith(prefix) and f.name.endswith("_dos.json"):
            full_prefix = f.name.rpartition("_")[0]
            with_prefix.append((str(f), full_prefix))

    return with_prefix

def _main():
    parser = argparse.ArgumentParser("Plot DOS(1) - DOS(2)")
    parser.add_argument("prefix", type=str, help="Calculation prefix")
    parser.add_argument("dir_add", type=str, help="Directory for DOS(1)")
    parser.add_argument("dir_sub", type=str, help="Directory for DOS(2)")
    args = parser.parse_args()

    paths_add = find_with_prefix(args.prefix, args.dir_add)
    paths_sub = find_with_prefix(args.prefix, args.dir_sub)

    for (dos_path_add, full_prefix_add), (dos_path_sub, full_prefix_sub) in zip(paths_add, paths_sub):
        assert(full_prefix_add == full_prefix_sub)

        with open(dos_path_add) as fp:
            dos_data_add = json.load(fp)

        with open(dos_path_sub) as fp:
            dos_data_sub = json.load(fp)

        es = dos_data_add["es"]
        total_dos_add = dos_data_add["total_dos"]
        total_dos_sub = dos_data_sub["total_dos"]
        ys = [add - sub for add, sub in zip(total_dos_add, total_dos_sub)]

        if "xlabel" in dos_data_add:
            plt.xlabel(dos_data_add["xlabel"])

        if "ylabel" in dos_data_add:
            plt.ylabel("$\\Delta$" + dos_data_add["ylabel"])

        if "caption" in dos_data_add:
            plt.title(dos_data_add["caption"] + " $- \\theta = {:.3}$ deg.".format(dos_data_sub["theta_deg"]))

        plot_path = "{}_dos_diff.png".format(full_prefix_add)

        plt.xlim(min(es), max(es))
        plt.ylim(min(ys), max(ys))

        plt.plot(es, ys, 'k-')

        plt.savefig(plot_path, dpi=500, bbox_inches='tight')
        plt.clf()

if __name__ == "__main__":
    _main()
