"""
Plot 2D maps produced by the `mesh` plan in scan_8id.py.

The `mesh` plan rasters motor1 (outer/slow) over motor2 (inner/fast) and stores
`shape: [num1, num2]` plus `plan_name: "mesh"` in the run start metadata. This
module reads a run back from the databroker catalog and displays the chosen
detector signal as a 2D image.

usage (in the bluesky session):
    from id8_common.plans.align.plot_mesh import plot_mesh
    plot_mesh()                                  # last scan, default signal
    plot_mesh(-1, signal="lambda2M_stats2_total")
    plot_mesh("b3a4aebf")                         # by (partial) uid or scan_id
"""

import numpy as np
import matplotlib.pyplot as plt

from apsbits.core.instrument_init import oregistry  # noqa: F401 (session context)

try:
    from databroker import catalog  # noqa: F401
except Exception:  # pragma: no cover - catalog is provided by the session
    pass


def _find_catalog():
    """Locate the session catalog `cat` from the IPython user namespace."""
    try:
        from IPython import get_ipython

        ip = get_ipython()
        if ip is not None:
            cat = ip.user_ns.get("cat")
            if cat is not None:
                return cat
    except Exception:
        pass
    return None


def _get_run(cat, ref):
    """Resolve a run from the catalog by index, uid, or scan_id."""
    if ref is None:
        ref = -1
    return cat[ref]


def plot_mesh(
    ref=-1,
    signal="lambda2M_stats2_total",
    cat=None,
    ax=None,
    cmap="viridis",
    log=False,
):
    """
    Plot a 2D map from a `mesh` run.

    Args:
        ref: run reference - integer index (-1 = last), uid string, or scan_id.
        signal: data column to image (e.g. "lambda2M_stats1_total",
            "lambda2M_stats2_total", "eiger4M_stats1_total", ...).
        cat: databroker catalog. If None, uses the session global `cat`.
        ax: existing matplotlib Axes to draw into (optional).
        cmap: matplotlib colormap name.
        log: if True, display log10 of the signal (clipped at >0).

    Returns:
        (fig, ax, image_2d) so the caller can further customize or re-save.
    """
    if cat is None:
        cat = _find_catalog()
        if cat is None:
            raise RuntimeError(
                "No catalog available. Pass cat=<catalog> explicitly."
            )

    run = _get_run(cat, ref)
    start = run.metadata["start"]

    if start.get("plan_name") != "mesh":
        print(f"WARNING: run plan_name is {start.get('plan_name')!r}, not 'mesh'.")

    ds = run.primary.read()

    if signal not in ds:
        avail = [k for k in ds.keys() if "total" in k or "stats" in k]
        raise KeyError(
            f"Signal {signal!r} not in run. Candidates: {avail}"
        )

    z = np.asarray(ds[signal]).ravel()

    # grid shape and axis names/positions
    shape = start.get("shape")
    motors = start.get("motors", ["motor1", "motor2"])
    m1, m2 = motors[0], motors[1]

    if shape is None:
        n = int(round(np.sqrt(z.size)))
        shape = [n, n]
        print(f"WARNING: no shape in metadata; assuming {shape}.")

    num1, num2 = int(shape[0]), int(shape[1])
    npts = num1 * num2

    # Place each measured point into its true grid cell by actual motor
    # position, rather than blindly reshaping the event stream. This is
    # robust to dropped/duplicated points, aborted scans, event ordering,
    # and non-square grids -- all of which corrupt a plain reshape.
    if m1 in ds and m2 in ds:
        x1 = np.asarray(ds[m1]).ravel()  # motor1 (outer/slow)
        x2 = np.asarray(ds[m2]).ravel()  # motor2 (inner/fast)
        n = min(len(x1), len(x2), len(z))
        x1, x2, zz = x1[:n], x2[:n], z[:n]

        # reconstruct the commanded grid axes from the measured span
        ax1 = np.linspace(np.nanmin(x1), np.nanmax(x1), num1)  # outer -> rows
        ax2 = np.linspace(np.nanmin(x2), np.nanmax(x2), num2)  # inner -> cols

        img = np.full((num1, num2), np.nan, dtype=float)
        for xi, yi, si in zip(x1, x2, zz):
            ir = int(np.argmin(np.abs(ax1 - xi)))
            ic = int(np.argmin(np.abs(ax2 - yi)))
            img[ir, ic] = si

        n_filled = int(np.isfinite(img).sum())
        if n_filled != npts:
            print(f"NOTE: {n_filled}/{npts} grid cells filled "
                  f"({z.size} events for a {num1}x{num2} map).")

        extent = [
            float(ax2.min()), float(ax2.max()),
            float(ax1.min()), float(ax1.max()),
        ]
    else:
        # fallback: no motor columns -> blind reshape with NaN padding
        if z.size < npts:
            z = np.concatenate([z[:npts], np.full(npts - z.size, np.nan)])
        else:
            z = z[:npts]
        img = z.reshape(num1, num2).astype(float)
        extent = None

    disp = img.astype(float)
    label = signal
    if log:
        disp = np.log10(np.clip(disp, 1, None))
        label = f"log10({signal})"

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.figure

    im = ax.imshow(
        disp,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        extent=extent,
        interpolation="nearest",
    )
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(label)

    ax.set_xlabel(f"{m2}")
    ax.set_ylabel(f"{m1}")
    scan_id = start.get("scan_id", "?")
    uid = start.get("uid", "")[:8]
    ax.set_title(f"mesh scan {scan_id} [{uid}]\n{signal}")

    fig.tight_layout()
    plt.show()
    return fig, ax, img
