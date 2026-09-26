#!/usr/bin/env python3
"""Live viewer for the LBLS solver's results/N dumps.

Reads the block-format text dumps (clArray3D::save2file(): nX blocks of
nY rows x nZ cols, blank-line separated; vector fields have nZ*3 cols
laid out as x0,y0,z0,x1,y1,z1,...) and renders a Y-midplane slice of
density plus in-plane velocity vectors, matching the combined view
LivePostProcessor.m used to produce.
"""
import sys
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_scalar_block(path, nX, nY, nZ):
    arr = np.zeros((nX, nY, nZ))
    with open(path) as f:
        blocks_text = f.read().split("\n\n")
    for xi, block in enumerate(blocks_text):
        block = block.strip("\n")
        if not block:
            continue
        rows = [line.split() for line in block.split("\n") if line.strip()]
        data = np.array(rows, dtype=float)
        arr[xi] = data[:nY, :nZ]
    return arr


def load_vector_block(path, nX, nY, nZ):
    arr = np.zeros((nX, nY, nZ, 3))
    with open(path) as f:
        blocks_text = f.read().split("\n\n")
    for xi, block in enumerate(blocks_text):
        block = block.strip("\n")
        if not block:
            continue
        rows = [line.split() for line in block.split("\n") if line.strip()]
        data = np.array(rows, dtype=float)
        data = data[:nY, : nZ * 3].reshape(nY, nZ, 3)
        arr[xi] = data
    return arr


def load_structure_positions(path):
    vals = []
    with open(path) as f:
        for line in f:
            vals.extend(float(x) for x in line.split())
    vals = np.array(vals).reshape(-1, 3)
    return vals


def latest_dump_dir(base="results"):
    dirs = [d for d in glob.glob(os.path.join(base, "*")) if os.path.isdir(d)]
    nums = []
    for d in dirs:
        b = os.path.basename(d)
        if b.isdigit():
            nums.append((int(b), d))
    if not nums:
        return None, None
    nums.sort()
    return nums[-1]


def main():
    with open("load/lbsize.txt") as f:
        nX, nY, nZ = (int(float(x)) for x in f.read().split())

    dump_num, dump_dir = latest_dump_dir()
    if dump_dir is None:
        print("No dumps yet.")
        sys.exit(1)

    print(f"Rendering dump {dump_num} from {dump_dir}")

    ro = load_scalar_block(os.path.join(dump_dir, "lbro.txt"), nX, nY, nZ)
    j = load_vector_block(os.path.join(dump_dir, "lbj.txt"), nX, nY, nZ)

    midY = nY // 2
    ro_slice = ro[:, midY, :].T  # shape (nZ, nX)
    jx_slice = j[:, midY, :, 0].T
    jz_slice = j[:, midY, :, 2].T

    fig, ax = plt.subplots(figsize=(12, 4))
    vmax = np.max(np.abs(ro_slice)) or 1e-12
    im = ax.imshow(
        ro_slice,
        origin="lower",
        extent=[0, nX, 0, nZ],
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label("density deviation")

    step_x = max(nX // 40, 1)
    step_z = max(nZ // 20, 1)
    Xg, Zg = np.meshgrid(
        np.arange(0, nX, step_x) + 0.5, np.arange(0, nZ, step_z) + 0.5
    )
    U = jx_slice[::step_z, ::step_x]
    W = jz_slice[::step_z, ::step_x]
    speed = np.sqrt(U**2 + W**2)
    if speed.max() > 0:
        ax.quiver(Xg, Zg, U, W, color="k", scale=speed.max() * 25, width=0.002)

    try:
        lsc = load_structure_positions(os.path.join(dump_dir, "lsc.txt"))
        lso = np.loadtxt("load/lso.txt")
        wall_mask = lso == 0
        wall = lsc[wall_mask]
        near_mid = np.abs(wall[:, 1] - midY) < 1.5
        ax.scatter(
            wall[near_mid, 0],
            wall[near_mid, 2],
            s=2,
            c="gray",
            alpha=0.5,
            label="vessel wall",
        )
        valve_mask = ~wall_mask
        valve = lsc[valve_mask]
        near_mid_v = np.abs(valve[:, 1] - midY) < 1.5
        ax.scatter(
            valve[near_mid_v, 0],
            valve[near_mid_v, 2],
            s=4,
            c="magenta",
            label="valve/leaflet",
        )
        ax.legend(loc="upper right", fontsize=8)
    except Exception as e:
        print("structure overlay skipped:", e)

    ax.set_xlabel("X (lattice units)")
    ax.set_ylabel("Z (lattice units)")
    ax.set_title(f"Lymphangion LBM/LSM — dump {dump_num} (Y-midplane slice)")

    out_path = "live_view_latest.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
