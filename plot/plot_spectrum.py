import argparse
import json
import matplotlib.pyplot as plt

def _main():
    parser = argparse.ArgumentParser("plot specta")
    parser.add_argument("prefix", type=str, help="Prefix for calculation")
    args = parser.parse_args()

    bands_path = "{}_spectrum.json".format(args.prefix)

    with open(bands_path, 'r') as fp:
        bands_data = json.load(fp)

    kxs = [k[0] for k in bands_data["ks"]]

    for em in bands_data["emks"]:
        plt.plot(kxs, em, 'k-')

    plt.show()

if __name__ == "__main__":
    _main()
