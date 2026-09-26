#!/usr/bin/env python3
"""Live multipanel viewer for the LBLS solver's results/N dumps.

Panel A (3D): vessel wall + valve/leaflet structure (colored by radial
deformation from t=0) with a coarse 3D fluid-velocity quiver overlay.
Panel B: inlet vs outlet mass flow rate over time (mass-flux balance
across the lymphangion, from the solver's own tout.txt diagnostic log).
Panel C: compressibility indicator over time (max / mean-abs density
deviation from the reference density — this solver is a weakly
compressible LBM, so density deviation is the direct compressibility
signal).
"""
import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def load_flat(path):
    vals = []
    with open(path) as f:
        for line in f:
            vals.extend(float(x) for x in line.split())
    return np.array(vals)


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


def latest_dump_dir(base="results"):
    dirs = [d for d in glob.glob(os.path.join(base, "*")) if os.path.isdir(d)]
    nums = [(int(os.path.basename(d)), d) for d in dirs if os.path.basename(d).isdigit()]
    if not nums:
        return None, None
    nums.sort()
    return nums[-1]


def load_tout(path="tout.txt"):
    """Columns: Time, mass, NM.S, Qavgl, Qavgm, Qavgh, avgl, avgm, avgh, Ro.max, Ro.avgabs, [LS...]"""
    rows = []
    if not os.path.exists(path):
        return None
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 11:
                continue
            rows.append([float(x) for x in parts[:11]])
    if not rows:
        return None
    return np.array(rows)


def main():
    nX, nY, nZ = (int(float(x)) for x in load_flat("load/lbsize.txt"))

    dump_num, dump_dir = latest_dump_dir()
    if dump_dir is None:
        print("No dumps yet.")
        sys.exit(1)
    print(f"Rendering dump {dump_num} from {dump_dir}")

    lso = load_flat("load/lso.txt")
    N = len(lso)
    c0 = load_flat("load/lsc.txt").reshape(N, 3)
    c_now = load_flat(os.path.join(dump_dir, "lsc.txt")).reshape(N, 3)

    midY, midZ = (nY - 1) / 2.0, (nZ - 1) / 2.0
    r0 = np.sqrt((c0[:, 1] - midY) ** 2 + (c0[:, 2] - midZ) ** 2)
    r_now = np.sqrt((c_now[:, 1] - midY) ** 2 + (c_now[:, 2] - midZ) ** 2)
    dr = r_now - r0

    j = load_vector_block(os.path.join(dump_dir, "lbj.txt"), nX, nY, nZ)

    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.4, 1], hspace=0.32, wspace=0.28)
    ax3d = fig.add_subplot(gs[:, 0], projection="3d")
    ax_flow = fig.add_subplot(gs[0, 1])
    ax_comp = fig.add_subplot(gs[1, 1])

    # --- Panel A: 3D structure + velocity ---
    wall = c_now[lso == 0]
    wall_ds = wall[:: max(len(wall) // 1500, 1)]
    ax3d.scatter(
        wall_ds[:, 0], wall_ds[:, 1], wall_ds[:, 2],
        s=3, c="gray", alpha=0.25, label="vessel wall",
    )

    valve_mask = lso != 0
    valve = c_now[valve_mask]
    valve_dr = dr[valve_mask]
    vmax = max(np.max(np.abs(valve_dr)), 1e-9)
    sc = ax3d.scatter(
        valve[:, 0], valve[:, 1], valve[:, 2],
        s=8, c=valve_dr, cmap="coolwarm", vmin=-vmax, vmax=vmax,
        label="valve/leaflet (color = radial deformation)",
    )
    cbar = fig.colorbar(sc, ax=ax3d, pad=0.08, shrink=0.6)
    cbar.set_label("leaflet Δr from t=0")

    sx, sy, sz = max(nX // 18, 1), max(nY // 6, 1), max(nZ // 6, 1)
    Xg, Yg, Zg = np.meshgrid(
        np.arange(0, nX, sx), np.arange(0, nY, sy), np.arange(0, nZ, sz), indexing="ij"
    )
    U = j[::sx, ::sy, ::sz, 0]
    V = j[::sx, ::sy, ::sz, 1]
    W = j[::sx, ::sy, ::sz, 2]
    speed = np.sqrt(U**2 + V**2 + W**2)
    if speed.max() > 0:
        ax3d.quiver(
            Xg, Yg, Zg, U, V, W,
            length=6.0 / speed.max(), normalize=False, color="k", alpha=0.6, linewidth=0.6,
        )

    ax3d.set_xlabel("X")
    ax3d.set_ylabel("Y")
    ax3d.set_zlabel("Z")
    ax3d.set_box_aspect((nX / nY, 1, 1))
    ax3d.set_title(f"Dump {dump_num} — structure deformation + flow field")
    ax3d.legend(loc="upper left", fontsize=7)

    # --- Panels B/C: inlet/outlet mass flow + compressibility over time ---
    tout = load_tout()
    if tout is not None:
        t = tout[:, 0]
        avgl, avgm, avgh = tout[:, 6], tout[:, 7], tout[:, 8]
        ro_max, ro_avgabs = tout[:, 9], tout[:, 10]

        ax_flow.plot(t, avgl, label="inlet (low X)", color="tab:blue")
        ax_flow.plot(t, avgm, label="mid X", color="gray", linestyle="--", linewidth=1)
        ax_flow.plot(t, avgh, label="outlet (high X)", color="tab:red")
        ax_flow.axhline(0, color="k", linewidth=0.5)
        ax_flow.set_ylabel("mass flow rate")
        ax_flow.set_title("Inlet vs outlet mass flow rate")
        ax_flow.legend(fontsize=8, loc="best")
        ax_flow.grid(alpha=0.3)

        ax_comp.plot(t, ro_max, label="max |Δρ| (peak compressibility)", color="tab:purple")
        ax_comp.plot(t, ro_avgabs, label="mean |Δρ|", color="tab:orange")
        ax_comp.set_xlabel("simulation time (lattice steps)")
        ax_comp.set_ylabel("density deviation")
        ax_comp.set_title("Compressibility (density deviation from ρ₀)")
        ax_comp.legend(fontsize=8, loc="best")
        ax_comp.grid(alpha=0.3)
    else:
        for ax, msg in [(ax_flow, "tout.txt not found"), (ax_comp, "tout.txt not found")]:
            ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes)

    fig.suptitle(
        f"Lymphangion LBM/LSM — dump {dump_num}", fontsize=13, fontweight="bold"
    )
    out_path = "live_view_latest.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
