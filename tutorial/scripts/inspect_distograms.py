#!/usr/bin/env python3
"""
inspect_distograms.py

Compute the average-distance matrix from a structure-predictor distogram
(.pkl: AF2 / ESMFold / other formats) and optionally overlay the residue pairs
restrained by a bAIes parameter file (.dat).

Examples
--------
    python inspect_distograms.py --pkl PaaA2_af2.pkl --out aves.npy
    python inspect_distograms.py --pkl PaaA2_af2.pkl --plot aves.pdf \
        --baies baies_PaaA2.dat --pdb PaaA2.pdb
"""
import argparse
import sys
from pickle import load

import numpy as np
from scipy.special import softmax, expit

# Standard AF2 bin edges (Angstrom), used when 'bin_edges' is missing (e.g. ESMFold2)
AF2_BIN_EDGES = np.array([
    2.3125, 2.625, 2.9375, 3.25, 3.5625, 3.875,
    4.1875, 4.5, 4.8125, 5.125, 5.4375, 5.75,
    6.0625, 6.375, 6.6875, 7., 7.3125, 7.625,
    7.9375, 8.25, 8.5625, 8.875, 9.1875, 9.5,
    9.812499, 10.125, 10.4375, 10.75, 11.0625, 11.375,
    11.6875, 12., 12.3125, 12.625, 12.9375, 13.25,
    13.5625, 13.875, 14.1875, 14.499999, 14.8125, 15.125,
    15.4375, 15.75, 16.0625, 16.375, 16.6875, 16.999998,
    17.312498, 17.625, 17.9375, 18.25, 18.5625, 18.875,
    19.1875, 19.5, 19.8125, 20.125, 20.437498, 20.75,
    21.0625, 21.375, 21.6875], dtype="float32")

# Value assigned to pairs whose probability mass is mostly in the >22 A bin (nm)
FAR_VALUE = 2.2
# Marker values for bAIes-restrained pairs (outside [vmin, vmax] -> drawn in red)
MARK_I4 = -1.0   # |i-j| == 4 pairs
MARK_OTHER = 2.3  # all other restrained pairs


def read_pkl(pkl_file, convf="softmax"):
    """
    Read a distogram .pkl file and return
      bins       : (n_bins-1,) bin centers in nm (last open-ended bin excluded)
      distograms : (N, N, n_bins) probability distributions (last bin = >22 A)
    """
    with open(pkl_file, "rb") as f:
        data = load(f)["distogram"]

    try:
        bins = np.asarray(data["bin_edges"], dtype=float)
    except KeyError:
        print("bin_edges not in the pkl file. Assuming standard AF2 format for pair rep...")
        bins = AF2_BIN_EDGES.astype(float)

    # Shift edges to bin centers and convert Angstrom -> nm
    delta_bins = bins[1] - bins[0]
    bins = (bins - delta_bins / 2) * 0.1

    logits = np.asarray(data["logits"])
    if logits.ndim == 4:
        # ESMFold .pkl format: (1, N, N, n_bins)
        logits = logits[0]

    if convf == "softmax":
        distograms = softmax(logits, axis=-1)
    elif convf == "sigmoid":
        distograms = expit(logits)
        distograms = distograms / distograms.sum(axis=-1, keepdims=True)
    else:
        sys.exit("Conversion function not recognized. Select softmax (default) or sigmoid.")

    if distograms.shape[-1] != len(bins) + 1:
        sys.exit("Inconsistent shapes: {} distogram bins but {} bin edges.".format(
            distograms.shape[-1], len(bins)))
    return bins, distograms


def average_distances(bins, distograms, far_threshold=0.1):
    """
    Expected distance (nm) per residue pair, computed over the finite bins only.
    Pairs with P(>22 A) >= far_threshold are set to FAR_VALUE.
    """
    finite = distograms[..., :-1]
    aves = (finite * bins).sum(axis=-1) / finite.sum(axis=-1)
    aves[distograms[..., -1] >= far_threshold] = FAR_VALUE
    return aves


def read_baies_pairs(baies_file, topology, n_header=2):
    """
    Read the restrained atom pairs from a bAIes .dat file (columns 2 and 3 are
    atom serials) and map them to 0-based residue indices using the topology.
    """
    try:
        import MDAnalysis as mda
    except ImportError:
        sys.exit("MDAnalysis is required to map bAIes atom indices to residues.")

    with open(baies_file, "r") as f:
        lines = [el.split() for el in f.readlines()][n_header:]
    lines = [el for el in lines if el]
    at1 = [int(el[1]) for el in lines]
    at2 = [int(el[2]) for el in lines]

    u = mda.Universe(topology)
    id2res = dict(zip(u.atoms.ids, u.atoms.resindices))
    try:
        pairs = {(int(id2res[a]), int(id2res[b])) for a, b in zip(at1, at2)}
    except KeyError as e:
        sys.exit("Atom id {} from {} not found in {}.".format(e, baies_file, topology))
    return pairs, u.residues.n_residues


def mark_pairs(aves, pairs):
    """Return a copy of aves with bAIes-restrained pairs set to marker values."""
    marked = aves.copy()
    for i, j in pairs:
        val = MARK_I4 if abs(i - j) == 4 else MARK_OTHER
        marked[i, j] = val
        marked[j, i] = val
    return marked


def plot_matrix(aves, out_pdf, title="Distograms average distance matrix"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(88 / 25.4, 70 / 25.4))

    cmap = matplotlib.colormaps.get_cmap("Blues").reversed()
    cmap.set_under("red")
    cmap.set_over("red")

    plotmap = ax.pcolormesh(aves, cmap=cmap, rasterized=True, vmin=0, vmax=FAR_VALUE)
    cbar = plt.colorbar(plotmap, ax=ax)
    cbar.set_label("Average distance (nm)")
    ax.set_title(title)
    ax.invert_yaxis()
    ax.set_xlabel("Residues")
    ax.set_ylabel("Residues")

    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot / save the average distance matrix from a distogram .pkl file.")
    parser.add_argument("--pkl", required=True, help="Distogram .pkl file.")
    parser.add_argument("--baies", help="bAIes parameter file (.dat). Restrained pairs are marked in red.")
    parser.add_argument("--pdb", help="Topology (.pdb) matching the bAIes atom indices. Required with --baies.")
    parser.add_argument("--out", help="Output .npy file for the average distance matrix (nm).")
    parser.add_argument("--plot", help="Output .pdf file for the plot.")
    parser.add_argument("--far-threshold", type=float, default=0.1,
                        help="If P(>22 A) >= this value, the pair is set to {} nm (default: 0.1).".format(FAR_VALUE))
    parser.add_argument("--title", default="Distograms average distance matrix", help="Plot title.")
    args = parser.parse_args()

    if args.out is None and args.plot is None:
        parser.error("at least one of --out (.npy) or --plot (.pdf) must be specified.")
    if args.baies and not args.pdb:
        parser.error("--pdb is required with --baies (to map atom indices to residues).")
    if args.pdb and not args.baies:
        print("Warning: --pdb given without --baies; it will be ignored.")
    return args


def main():
    args = parse_args()

    bins, distograms = read_pkl(args.pkl, convf='softmax')
    aves = average_distances(bins, distograms, far_threshold=args.far_threshold)
    print("Average distance matrix: {} x {} residues".format(*aves.shape))

    if args.out:
        np.save(args.out, aves)
        print("Saved matrix to {}".format(args.out))

    if args.plot:
        to_plot = aves
        if args.baies:
            pairs, n_res = read_baies_pairs(args.baies, args.pdb)
            if n_res != aves.shape[0]:
                print("Warning: topology has {} residues but distogram has {}.".format(n_res, aves.shape[0]))
            out_of_range = [p for p in pairs if max(p) >= aves.shape[0]]
            if out_of_range:
                sys.exit("{} bAIes pairs fall outside the distogram matrix.".format(len(out_of_range)))
            to_plot = mark_pairs(aves, pairs)
            print("Marked {} bAIes restrained pairs".format(len(pairs)))
        plot_matrix(to_plot, args.plot, title=args.title)
        print("Saved plot to {}".format(args.plot))


if __name__ == "__main__":
    main()
