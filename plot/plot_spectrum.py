import json
import matplotlib.pyplot as plt

def _main():
    prefix = "blg_moire"

    bands_path = "{}_spectrum.json".format(prefix)

    with open(bands_path, 'r') as fp:
        bands_data = json.load(fp)

    kxs = [k[0] for k in bands_data["ks"]]

    for em in bands_data["emks"]:
        plt.plot(kxs, em, 'k-')

    plt.show()

if __name__ == "__main__":
    _main()
