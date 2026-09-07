#!/usr/bin/env python
"""Render one 2D slice of /xpcs/multitau/unnormalized_G2 as a standalone
interactive HTML page (canvas + adjustable color range, no external JS).

Example
-------
    python scripts/plot_unnormalized_g2.py \
        /home/8-id-e/2026-2/pope202607/analysis/Multitau/A0096_..._G2_results.hdf \
        --tau -1 --qbin 0 -o /tmp/G2_slice.html
"""

import argparse
import base64
import json
from pathlib import Path

import h5py
import numpy as np
from matplotlib import colormaps
from scipy.ndimage import median_filter

DATASET = "/xpcs/multitau/unnormalized_G2"
CMAPS = ["viridis", "cividis", "inferno", "magma", "gray", "turbo"]


def build_colormaps(n=256):
    out = {}
    for name in CMAPS:
        rgba = colormaps[name](np.linspace(0, 1, n))
        out[name] = (rgba[:, :3] * 255).round().astype(np.uint8).ravel().tolist()
    return out


def find_hot_pixels(saxs, mask, ratio_thresh=20.0, floor=0.05, box=7):
    """Pixels far above their local neighbourhood in the time-averaged image."""
    local = median_filter(saxs, size=box)
    ratio = saxs / np.maximum(local, 1e-9)
    hot = mask & (saxs > floor) & (ratio > ratio_thresh)
    return np.argwhere(hot), local, ratio


def g2_by_static_bin(G2, IP, IF, sq_flat, nsq, counts):
    """boost_corr reduction: average pixels inside a static bin, then divide."""
    sums = np.stack([np.bincount(sq_flat, weights=a, minlength=nsq + 1) for a in (G2, IP, IF)])
    with np.errstate(invalid="ignore", divide="ignore"):
        avg = sums / counts
    ipif = avg[1] * avg[2]
    ipif = np.where(ipif == 0, 1e8, ipif)
    return avg[0] / ipif, sums


def hot_pixel_table(f, G2, IP, IF, hot_yx):
    """Per-pixel leverage on its dynamic-ROI g2, using the exact pipeline reduction."""
    dq2 = f["xpcs/qmap/dynamic_roi_map"][()]
    sq2 = f["xpcs/qmap/static_roi_map"][()]
    saxs = f["xpcs/temporal_mean/scattering_2d"][0].astype(np.float64)
    qlist = f["xpcs/qmap/dynamic_v_list_dim0"][()]
    dq, sq = dq2.ravel(), sq2.ravel()
    nsq = int(sq.max())
    counts = np.bincount(sq, minlength=nsq + 1).astype(np.float64)
    g2f, ipf, iff = G2.ravel(), IP.ravel(), IF.ravel()
    g2s, sums = g2_by_static_bin(g2f, ipf, iff, sq, nsq, counts)
    subbins = {i: np.unique(sq[dq == i]) for i in range(1, int(dq.max()) + 1)}
    subbins = {i: v[v > 0] for i, v in subbins.items()}
    g2_roi = {i: float(np.nanmean(g2s[v])) for i, v in subbins.items()}

    nx = dq2.shape[1]
    rows = []
    for y, x in hot_yx:
        j = y * nx + x
        s, dqi = int(sq[j]), int(dq[j])
        n = counts[s]
        new = np.nan
        if n > 1:
            a = (sums[0][s] - g2f[j]) / (n - 1)
            b = (sums[1][s] - ipf[j]) / (n - 1)
            c = (sums[2][s] - iff[j]) / (n - 1)
            new = a / (b * c) if b * c != 0 else np.nan
        nsub = len(subbins.get(dqi, [])) or 1
        rows.append({
            "y": int(y), "x": int(x), "dq": dqi, "sq": s,
            "I": float(saxs[y, x]),
            "IP": float(ipf[j]), "IF": float(iff[j]), "G2": float(g2f[j]),
            "g2s": float(g2s[s]), "g2s_wo": float(new),
            "lev": float((g2s[s] - new) / nsub) if np.isfinite(new) else 0.0,
        })
    rows.sort(key=lambda r: -r["lev"])

    # ROI g2 with every flagged pixel removed
    keep = np.ones(g2f.size, bool)
    keep[[r["y"] * nx + r["x"] for r in rows]] = False
    g2s_c, _ = g2_by_static_bin(
        g2f[keep], ipf[keep], iff[keep], sq[keep], nsq,
        np.bincount(sq[keep], minlength=nsq + 1).astype(np.float64))
    rois = [{
        "dq": i,
        "q": float(qlist[i - 1]),
        "npix": int((dq == i).sum()),
        "nhot": sum(1 for r in rows if r["dq"] == i),
        "g2": g2_roi[i],
        "g2c": float(np.nanmean(g2s_c[v])),
    } for i, v in subbins.items()]
    return rows, rois, dq2


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("hdf", type=Path)
    p.add_argument("--dataset", default=DATASET)
    p.add_argument("--tau", type=int, default=-1, help="index along axis 0 (default: last)")
    p.add_argument("--qbin", type=int, default=0, help="index along axis 1 (default: first)")
    p.add_argument("-o", "--output", type=Path, default=None)
    p.add_argument("--no-flag", action="store_true", help="skip hot-pixel detection")
    p.add_argument("--ratio", type=float, default=20.0, help="hot-pixel local-contrast threshold")
    args = p.parse_args()

    hot_rows, roi_rows, dq_map = [], [], None
    with h5py.File(args.hdf, "r") as f:
        dset = f[args.dataset]
        shape = dset.shape
        i0 = args.tau % shape[0]
        i1 = args.qbin % shape[1]
        img = np.asarray(dset[i0, i1, :, :], dtype=np.float32)
        delays = f["xpcs/multitau/delay_list"][()]
        if not args.no_flag:
            saxs = f["xpcs/temporal_mean/scattering_2d"][0].astype(np.float64)
            mask = f["xpcs/qmap/mask"][()]
            hot_yx, _, _ = find_hot_pixels(saxs, mask, ratio_thresh=args.ratio)
            G2, IP, IF = (dset[i0, k].astype(np.float64) for k in range(3))
            hot_rows, roi_rows, dq_map = hot_pixel_table(f, G2, IP, IF, hot_yx)
            del G2, IP, IF

    ny, nx = img.shape
    finite = np.isfinite(img)
    vals = img[finite]
    vmin_data = float(vals.min())
    vmax_data = float(vals.max())
    pos = vals[vals > 0]
    pos_min = float(pos.min()) if pos.size else 1e-6

    pcts = [0.0, 0.1, 1.0, 5.0, 50.0, 95.0, 99.0, 99.9, 100.0]
    pct_all = dict(zip(map(str, pcts), map(float, np.percentile(vals, pcts))))
    pct_pos = dict(zip(map(str, pcts), map(float, np.percentile(pos, pcts)))) if pos.size else pct_all

    # log-spaced histogram of positive values, for the range picker
    nbins = 256
    edges = np.logspace(np.log10(pos_min), np.log10(max(vmax_data, pos_min * 10)), nbins + 1)
    counts, _ = np.histogram(pos, bins=edges)
    n_zero = int((vals == 0).sum())
    n_nan = int((~finite).sum())

    payload = {
        "file": str(args.hdf),
        "dataset": args.dataset,
        "shape": list(shape),
        "i0": i0,
        "i1": i1,
        "nx": nx,
        "ny": ny,
        "vminData": vmin_data,
        "vmaxData": vmax_data,
        "posMin": pos_min,
        "pctAll": pct_all,
        "pctPos": pct_pos,
        "histEdges": edges.tolist(),
        "histCounts": counts.tolist(),
        "nZero": n_zero,
        "nNaN": n_nan,
        "nPix": int(img.size),
        "cmaps": build_colormaps(),
        "plane": ["G2", "IP", "IF"][i1] if i1 < 3 else str(i1),
        "delay": int(delays[i0]),
        "hot": hot_rows,
        "rois": roi_rows,
    }

    b64 = base64.b64encode(np.ascontiguousarray(img, dtype="<f4").tobytes()).decode("ascii")
    dq64 = ""
    if dq_map is not None:
        dq64 = base64.b64encode(np.ascontiguousarray(dq_map, dtype=np.uint8).tobytes()).decode("ascii")

    out = args.output or args.hdf.with_suffix("").with_name(
        f"{args.hdf.stem}_unnormG2_t{i0}_q{i1}.html"
    )
    out.write_text(
        HTML.replace("__META__", json.dumps(payload))
        .replace("__DATA__", b64)
        .replace("__DQMAP__", dq64)
    )
    print(f"slice [{i0}, {i1}] ({payload['plane']}, delay {payload['delay']}) -> {ny} x {nx}")
    print(f"value range: {vmin_data:g} .. {vmax_data:g}  (zeros: {n_zero}, NaN: {n_nan})")
    if hot_rows:
        print(f"flagged {len(hot_rows)} hot pixels; top leverage:")
        for r in hot_rows[:6]:
            print(f"   (row {r['y']:4d}, col {r['x']:4d})  ROI {r['dq']:2d}  "
                  f"I={r['I']:.4f} ct/frame  dg2_ROI={r['lev']:+.4f}")
    print(f"wrote {out}  ({out.stat().st_size / 1e6:.1f} MB)")


HTML = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>unnormalized_G2 slice</title>
<style>
:root{
  --surface:#ffffff; --surface-2:#f5f6f7; --ink:#1a1c1e; --ink-2:#4b5157; --ink-3:#7a8188;
  --line:#dfe2e5; --accent:#3d5afe;
}
@media (prefers-color-scheme: dark){
  :root{ --surface:#141618; --surface-2:#1c1f22; --ink:#e9ecef; --ink-2:#b3b9bf;
         --ink-3:#7f868d; --line:#2c3033; --accent:#8c9eff; }
}
*{box-sizing:border-box}
body{margin:0;background:var(--surface);color:var(--ink);
  font:13px/1.45 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
header{padding:14px 18px 10px;border-bottom:1px solid var(--line)}
h1{margin:0 0 3px;font-size:15px;font-weight:600;letter-spacing:-.01em}
.sub{color:var(--ink-3);font-size:11.5px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
  word-break:break-all}
main{display:flex;gap:16px;padding:14px 18px;align-items:flex-start;flex-wrap:wrap}
#stage{position:relative;flex:1 1 640px;min-width:340px;height:74vh;min-height:420px;
  background:var(--surface-2);border:1px solid var(--line);border-radius:6px;overflow:hidden;
  cursor:grab}
#stage.drag{cursor:grabbing}
#view{position:absolute;top:0;left:0;image-rendering:pixelated}
#hud{position:absolute;left:8px;bottom:8px;padding:4px 8px;border-radius:4px;
  background:color-mix(in srgb,var(--surface) 88%,transparent);border:1px solid var(--line);
  font:11.5px ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--ink-2);
  pointer-events:none;white-space:pre}
aside{flex:0 0 320px;display:flex;flex-direction:column;gap:14px}
.panel{border:1px solid var(--line);border-radius:6px;padding:12px}
.panel h2{margin:0 0 9px;font-size:11px;font-weight:600;letter-spacing:.06em;
  text-transform:uppercase;color:var(--ink-3)}
label{display:block;font-size:11.5px;color:var(--ink-2);margin:9px 0 3px}
label:first-of-type{margin-top:0}
select,input[type=number]{width:100%;padding:5px 7px;font-size:12.5px;color:var(--ink);
  background:var(--surface);border:1px solid var(--line);border-radius:4px;
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
input[type=range]{width:100%;margin:2px 0;accent-color:var(--accent)}
.row{display:flex;gap:8px}.row>*{flex:1}
.btns{display:flex;flex-wrap:wrap;gap:6px;margin-top:10px}
button{padding:4px 9px;font-size:11.5px;color:var(--ink-2);background:var(--surface-2);
  border:1px solid var(--line);border-radius:4px;cursor:pointer}
button:hover{color:var(--ink);border-color:var(--ink-3)}
#hist{display:block;width:100%;height:64px;margin-top:4px;border-radius:3px;background:var(--surface-2)}
.cbwrap{display:flex;gap:8px;align-items:stretch;margin-top:2px}
#cbar{width:26px;height:132px;border:1px solid var(--line);border-radius:3px}
.cbticks{display:flex;flex-direction:column;justify-content:space-between;
  font:11px ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--ink-3)}
.stats{font:11.5px/1.6 ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--ink-2)}
.stats b{font-weight:600;color:var(--ink)}
.note{font-size:11px;color:var(--ink-3);margin-top:8px}
#overlay{position:absolute;top:0;left:0;pointer-events:none}
table.hot{width:100%;border-collapse:collapse;
  font:11px ui-monospace,SFMono-Regular,Menlo,monospace}
table.hot th{text-align:right;font-weight:600;color:var(--ink-3);padding:2px 4px;
  border-bottom:1px solid var(--line);white-space:nowrap}
table.hot td{text-align:right;padding:2px 4px;color:var(--ink-2);white-space:nowrap}
table.hot tbody tr{cursor:pointer}
table.hot tbody tr:hover td{background:var(--surface-2);color:var(--ink)}
table.hot tbody tr.sel td{background:color-mix(in srgb,var(--accent) 16%,transparent);color:var(--ink)}
.scroll{max-height:230px;overflow:auto}
.tag{display:inline-block;padding:1px 5px;border-radius:3px;font-size:10.5px;
  background:var(--surface-2);border:1px solid var(--line);color:var(--ink-2)}
</style></head><body>
<header>
  <h1>unnormalized_G2 &mdash; <span id="ttl"></span></h1>
  <div class="sub" id="src"></div>
</header>
<main>
  <div id="stage"><canvas id="view"></canvas><canvas id="overlay"></canvas>
    <div id="hud">move cursor over image</div></div>
  <aside>
    <div class="panel">
      <h2>Color range</h2>
      <label for="smin">vmin <span id="lmin" style="float:right;font-family:ui-monospace"></span></label>
      <input type="range" id="smin" min="0" max="10000" value="0">
      <label for="smax">vmax <span id="lmax" style="float:right;font-family:ui-monospace"></span></label>
      <input type="range" id="smax" min="0" max="10000" value="10000">
      <div class="row" style="margin-top:8px">
        <div><label for="nmin">min</label><input type="number" id="nmin" step="any"></div>
        <div><label for="nmax">max</label><input type="number" id="nmax" step="any"></div>
      </div>
      <canvas id="hist" width="600" height="128"></canvas>
      <div class="btns">
        <button data-p="0,100">full</button>
        <button data-p="1,99">1&ndash;99%</button>
        <button data-p="5,95">5&ndash;95%</button>
        <button data-p="0.1,99.9">0.1&ndash;99.9%</button>
        <button data-p="50,99.9">50&ndash;99.9%</button>
      </div>
      <div class="note">Sliders are log-spaced: fine control near zero. Percentile
        buttons use non-zero pixels only.</div>
    </div>
    <div class="panel">
      <h2>Mapping</h2>
      <div class="row">
        <div><label for="cmap">colormap</label>
          <select id="cmap"></select></div>
        <div><label for="scale">scale</label>
          <select id="scale">
            <option value="linear">linear</option>
            <option value="log" selected>log</option>
            <option value="sqrt">sqrt</option>
          </select></div>
      </div>
      <label style="margin-top:9px">colorbar</label>
      <div class="cbwrap">
        <canvas id="cbar" width="26" height="132"></canvas>
        <div class="cbticks"><span id="cbhi"></span><span id="cbmid"></span><span id="cblo"></span></div>
      </div>
      <div class="btns">
        <button id="rev">reverse map</button>
        <button id="flipy">flip Y</button>
        <button id="fit">fit view</button>
        <button id="png">save PNG</button>
      </div>
    </div>
    <div class="panel" id="hotpanel">
      <h2>Flagged hot pixels</h2>
      <div class="btns" style="margin:0 0 8px">
        <button id="tflag">markers: on</button>
        <button id="tlabel">labels: off</button>
        <button id="prev">&larr; prev</button>
        <button id="next">next &rarr;</button>
      </div>
      <div class="scroll"><table class="hot"><thead><tr>
        <th>row</th><th>col</th><th>ROI</th><th>ct/frm</th><th>&Delta;g2</th></tr></thead>
        <tbody id="hotbody"></tbody></table></div>
      <div class="note" id="hotnote"></div>
    </div>
    <div class="panel" id="roipanel">
      <h2>ROI g2 at this delay</h2>
      <div class="scroll"><table class="hot"><thead><tr>
        <th>ROI</th><th>q (&Aring;<sup>-1</sup>)</th><th>hot</th><th>g2</th><th>g2 clean</th></tr></thead>
        <tbody id="roibody"></tbody></table></div>
      <div class="note">"g2 clean" = the same reduction with every flagged pixel
        dropped from its static bin.</div>
    </div>
    <div class="panel">
      <h2>Slice</h2>
      <div class="stats" id="stats"></div>
    </div>
  </aside>
</main>
<script id="meta" type="application/json">__META__</script>
<script id="blob" type="text/plain">__DATA__</script>
<script id="dqmap" type="text/plain">__DQMAP__</script>
<script>
const M = JSON.parse(document.getElementById('meta').textContent);
const raw = atob(document.getElementById('blob').textContent.trim());
const bytes = new Uint8Array(raw.length);
for (let i = 0; i < raw.length; i++) bytes[i] = raw.charCodeAt(i);
const IMG = new Float32Array(bytes.buffer);
const NX = M.nx, NY = M.ny;

const dqRaw = document.getElementById('dqmap').textContent.trim();
let DQ = null;
if (dqRaw){
  const s = atob(dqRaw); DQ = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) DQ[i] = s.charCodeAt(i);
}
const HOT = M.hot || [];

document.getElementById('ttl').textContent =
  `[${M.i0}, ${M.i1}, :, :] = ${M.plane}, delay ${M.delay}  (${NY} × ${NX})`;
document.getElementById('src').textContent = M.file + '  →  ' + M.dataset +
  '  ' + JSON.stringify(M.shape);
document.getElementById('stats').innerHTML =
  `axis0 index &nbsp;<b>${M.i0}</b> of ${M.shape[0]}<br>` +
  `axis1 index &nbsp;<b>${M.i1}</b> of ${M.shape[1]}<br>` +
  `data min &nbsp;&nbsp;&nbsp;<b>${fmt(M.vminData)}</b><br>` +
  `data max &nbsp;&nbsp;&nbsp;<b>${fmt(M.vmaxData)}</b><br>` +
  `min &gt; 0 &nbsp;&nbsp;&nbsp;&nbsp;<b>${fmt(M.posMin)}</b><br>` +
  `median&ne;0 &nbsp;&nbsp;<b>${fmt(M.pctPos['50.0'])}</b><br>` +
  `zeros &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>${(100*M.nZero/M.nPix).toFixed(1)}%</b>` +
  `&nbsp;&nbsp;NaN <b>${M.nNaN}</b>`;

function fmt(v){
  if (!isFinite(v)) return String(v);
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) return v.toExponential(3);
  return v.toPrecision(a >= 1 ? 6 : 4).replace(/\.?0+$/, '');
}

/* ---- log-spaced slider mapping over [0, vmaxData] ---------------------- */
const K = 9, SMAX = 10000;
const pos2val = p => M.vmaxData * Math.expm1(K * p / SMAX) / Math.expm1(K);
const val2pos = v => Math.round(SMAX / K * Math.log1p(Math.max(0, v) / M.vmaxData * Math.expm1(K)));

/* ---- colormaps --------------------------------------------------------- */
const cmapSel = document.getElementById('cmap');
for (const k of Object.keys(M.cmaps)) {
  const o = document.createElement('option'); o.value = o.textContent = k; cmapSel.appendChild(o);
}
cmapSel.value = 'viridis';
let reversed = false, flipY = false;
function lut(){
  const a = M.cmaps[cmapSel.value], n = a.length / 3, out = new Uint8Array(a.length);
  for (let i = 0; i < n; i++){
    const j = reversed ? (n - 1 - i) : i;
    out[3*i] = a[3*j]; out[3*i+1] = a[3*j+1]; out[3*i+2] = a[3*j+2];
  }
  return out;
}

/* ---- state ------------------------------------------------------------- */
let vmin = M.pctPos['1.0'], vmax = M.pctPos['99.9'];
if (!(vmax > vmin)) { vmin = M.vminData; vmax = M.vmaxData; }

/* ---- offscreen render at native resolution ----------------------------- */
const off = document.createElement('canvas'); off.width = NX; off.height = NY;
const offCtx = off.getContext('2d');
const imgData = offCtx.createImageData(NX, NY);

function renderImage(){
  const cm = lut(), NC = cm.length / 3 - 1, d = imgData.data;
  const mode = document.getElementById('scale').value;
  const lo = (mode === 'log') ? Math.log(Math.max(vmin, M.posMin * 0.5)) : vmin;
  const hi = (mode === 'log') ? Math.log(Math.max(vmax, Math.max(vmin, M.posMin) * 1.000001)) : vmax;
  const span = (hi - lo) || 1;
  const floorV = Math.max(vmin, M.posMin * 0.5);
  for (let y = 0; y < NY; y++){
    const sy = flipY ? (NY - 1 - y) : y;
    let si = sy * NX, di = y * NX * 4;
    for (let x = 0; x < NX; x++, si++, di += 4){
      const v = IMG[si];
      let t;
      if (v !== v) { d[di] = d[di+1] = d[di+2] = 0; d[di+3] = 0; continue; }
      if (mode === 'log')       t = (Math.log(Math.max(v, floorV)) - lo) / span;
      else if (mode === 'sqrt') t = Math.sqrt(Math.max(0, (v - lo) / span));
      else                      t = (v - lo) / span;
      t = t < 0 ? 0 : (t > 1 ? 1 : t);
      const c = 3 * ((t * NC) | 0);
      d[di] = cm[c]; d[di+1] = cm[c+1]; d[di+2] = cm[c+2]; d[di+3] = 255;
    }
  }
  offCtx.putImageData(imgData, 0, 0);
  blit();
}

/* ---- view transform (zoom / pan) --------------------------------------- */
const stage = document.getElementById('stage');
const view = document.getElementById('view');
const ovl = document.getElementById('overlay');
const vctx = view.getContext('2d');
const octx = ovl.getContext('2d');
let z = 1, tx = 0, ty = 0;
let showFlags = true, showLabels = false, selHot = -1;

function resize(){
  const r = stage.getBoundingClientRect(), dpr = devicePixelRatio || 1;
  for (const c of [view, ovl]){
    c.width = Math.round(r.width * dpr); c.height = Math.round(r.height * dpr);
    c.style.width = r.width + 'px'; c.style.height = r.height + 'px';
    c.getContext('2d').setTransform(dpr, 0, 0, dpr, 0, 0);
  }
  blit();
}
function fitView(){
  const r = stage.getBoundingClientRect();
  z = Math.min(r.width / NX, r.height / NY) * 0.97;
  tx = (r.width - NX * z) / 2; ty = (r.height - NY * z) / 2;
  blit();
}
function blit(){
  const r = stage.getBoundingClientRect();
  vctx.clearRect(0, 0, r.width, r.height);
  vctx.imageSmoothingEnabled = z < 1;
  vctx.drawImage(off, tx, ty, NX * z, NY * z);
  drawFlags();
}
function drawFlags(){
  const r = stage.getBoundingClientRect();
  octx.clearRect(0, 0, r.width, r.height);
  if (!showFlags || !HOT.length) return;
  octx.font = '11px ui-monospace,SFMono-Regular,Menlo,monospace';
  octx.textBaseline = 'middle';
  HOT.forEach((h, i) => {
    const yy = flipY ? (NY - 1 - h.y) : h.y;
    const cx = tx + (h.x + 0.5) * z, cy = ty + (yy + 0.5) * z;
    if (cx < -40 || cy < -40 || cx > r.width + 40 || cy > r.height + 40) return;
    const sel = i === selHot;
    const rad = Math.max(7, Math.min(26, z * 3.5));
    octx.strokeStyle = sel ? '#ffd54f' : '#ff5252';
    octx.lineWidth = sel ? 2.5 : 1.6;
    octx.beginPath(); octx.arc(cx, cy, rad, 0, 6.2832); octx.stroke();
    octx.beginPath();
    octx.moveTo(cx - rad - 5, cy); octx.lineTo(cx - rad - 1, cy);
    octx.moveTo(cx + rad + 1, cy); octx.lineTo(cx + rad + 5, cy);
    octx.stroke();
    if (showLabels || sel){
      const t = `${h.y},${h.x}  ROI ${h.dq}  Δg2 ${h.lev >= 0 ? '+' : ''}${h.lev.toFixed(3)}`;
      const w = octx.measureText(t).width + 8;
      octx.fillStyle = 'rgba(0,0,0,.72)';
      octx.fillRect(cx + rad + 7, cy - 9, w, 18);
      octx.fillStyle = sel ? '#ffd54f' : '#ffb3b3';
      octx.fillText(t, cx + rad + 11, cy + 1);
    }
  });
}
stage.addEventListener('wheel', e => {
  e.preventDefault();
  const r = stage.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
  const f = Math.exp(-e.deltaY * 0.0016), nz = Math.min(60, Math.max(0.05, z * f));
  tx = mx - (mx - tx) * (nz / z); ty = my - (my - ty) * (nz / z); z = nz; blit();
}, {passive:false});
let drag = null;
stage.addEventListener('pointerdown', e => {
  drag = {x:e.clientX, y:e.clientY, tx, ty}; stage.classList.add('drag');
  stage.setPointerCapture(e.pointerId);
});
stage.addEventListener('pointerup', e => {
  drag = null; stage.classList.remove('drag'); stage.releasePointerCapture(e.pointerId);
});
stage.addEventListener('pointermove', e => {
  const r = stage.getBoundingClientRect();
  if (drag){ tx = drag.tx + (e.clientX - drag.x); ty = drag.ty + (e.clientY - drag.y); blit(); }
  const px = Math.floor(((e.clientX - r.left) - tx) / z);
  const py = Math.floor(((e.clientY - r.top) - ty) / z);
  const hud = document.getElementById('hud');
  if (px >= 0 && px < NX && py >= 0 && py < NY){
    const row = flipY ? (NY - 1 - py) : py;
    const k = row * NX + px;
    const roi = DQ ? DQ[k] : null;
    hud.textContent = `row ${row}  col ${px}   ${M.plane} ${fmt(IMG[k])}` +
      (roi !== null ? `   ROI ${roi || '—'}` : '') + `   zoom ${z.toFixed(2)}×`;
  } else hud.textContent = `zoom ${z.toFixed(2)}×`;
});

/* ---- colorbar + histogram ---------------------------------------------- */
function drawCbar(){
  const c = document.getElementById('cbar'), g = c.getContext('2d'), cm = lut(), NC = cm.length/3 - 1;
  const im = g.createImageData(c.width, c.height);
  for (let y = 0; y < c.height; y++){
    const t = 1 - y / (c.height - 1), j = 3 * ((t * NC) | 0);
    for (let x = 0; x < c.width; x++){
      const i = 4 * (y * c.width + x);
      im.data[i] = cm[j]; im.data[i+1] = cm[j+1]; im.data[i+2] = cm[j+2]; im.data[i+3] = 255;
    }
  }
  g.putImageData(im, 0, 0);
  const mode = document.getElementById('scale').value;
  let mid;
  if (mode === 'log'){
    const a = Math.log(Math.max(vmin, M.posMin*0.5)), b = Math.log(Math.max(vmax, M.posMin));
    mid = Math.exp((a + b) / 2);
  } else if (mode === 'sqrt'){ mid = vmin + (vmax - vmin) * 0.25; }
  else { mid = (vmin + vmax) / 2; }
  document.getElementById('cbhi').textContent = fmt(vmax);
  document.getElementById('cbmid').textContent = fmt(mid);
  document.getElementById('cblo').textContent = fmt(vmin);
}
function drawHist(){
  const c = document.getElementById('hist'), g = c.getContext('2d');
  const W = c.width, H = c.height, E = M.histEdges, C = M.histCounts;
  const cs = getComputedStyle(document.body);
  g.clearRect(0, 0, W, H);
  const lx0 = Math.log10(E[0]), lx1 = Math.log10(E[E.length-1]), lspan = lx1 - lx0;
  const cmax = Math.max(...C), lcmax = Math.log10(cmax + 1);
  g.fillStyle = cs.getPropertyValue('--ink-3').trim() || '#888';
  for (let i = 0; i < C.length; i++){
    if (!C[i]) continue;
    const x0 = W * (Math.log10(E[i]) - lx0) / lspan;
    const x1 = W * (Math.log10(E[i+1]) - lx0) / lspan;
    const h = H * Math.log10(C[i] + 1) / lcmax;
    g.fillRect(x0, H - h, Math.max(1, x1 - x0), h);
  }
  g.strokeStyle = cs.getPropertyValue('--accent').trim() || '#3d5afe';
  g.lineWidth = 2;
  for (const v of [vmin, vmax]){
    const x = W * (Math.log10(Math.max(v, E[0])) - lx0) / lspan;
    g.beginPath(); g.moveTo(x, 0); g.lineTo(x, H); g.stroke();
  }
}

/* ---- controls ---------------------------------------------------------- */
const smin = document.getElementById('smin'), smax = document.getElementById('smax');
const nmin = document.getElementById('nmin'), nmax = document.getElementById('nmax');
let pending = false;
function refresh(fromSlider){
  if (vmax <= vmin) vmax = vmin + Math.max(1e-12, Math.abs(vmin) * 1e-6);
  if (!fromSlider){ smin.value = val2pos(vmin); smax.value = val2pos(vmax); }
  nmin.value = Number(vmin.toPrecision(8)); nmax.value = Number(vmax.toPrecision(8));
  document.getElementById('lmin').textContent = fmt(vmin);
  document.getElementById('lmax').textContent = fmt(vmax);
  drawCbar(); drawHist();
  if (!pending){ pending = true; requestAnimationFrame(() => { pending = false; renderImage(); }); }
}
smin.addEventListener('input', () => {
  vmin = pos2val(+smin.value);
  if (vmin >= vmax){ vmax = pos2val(Math.min(10000, +smin.value + 1)); smax.value = val2pos(vmax); }
  refresh(true);
});
smax.addEventListener('input', () => {
  vmax = pos2val(+smax.value);
  if (vmax <= vmin){ vmin = pos2val(Math.max(0, +smax.value - 1)); smin.value = val2pos(vmin); }
  refresh(true);
});
nmin.addEventListener('change', () => { vmin = +nmin.value; refresh(false); });
nmax.addEventListener('change', () => { vmax = +nmax.value; refresh(false); });
for (const b of document.querySelectorAll('[data-p]')){
  b.addEventListener('click', () => {
    const [a, z2] = b.dataset.p.split(',');
    const src = (a === '0' && z2 === '100') ? M.pctAll : M.pctPos;
    vmin = src[(+a).toFixed(1)]; vmax = src[(+z2).toFixed(1)];
    if (a === '0' && z2 === '100'){ vmin = M.vminData; vmax = M.vmaxData; }
    refresh(false);
  });
}
/* ---- hot-pixel table --------------------------------------------------- */
function buildHot(){
  const tb = document.getElementById('hotbody');
  if (!HOT.length){
    document.getElementById('hotpanel').style.display = 'none';
    document.getElementById('roipanel').style.display = 'none';
    return;
  }
  HOT.forEach((h, i) => {
    const tr = document.createElement('tr');
    tr.innerHTML = `<td>${h.y}</td><td>${h.x}</td><td>${h.dq}</td>` +
      `<td>${h.I.toFixed(3)}</td><td>${h.lev >= 0 ? '+' : ''}${h.lev.toFixed(4)}</td>`;
    tr.addEventListener('click', () => gotoHot(i));
    tb.appendChild(tr);
  });
  const tot = HOT.reduce((a, h) => a + h.lev, 0);
  document.getElementById('hotnote').textContent =
    `${HOT.length} pixels > 20× their local 7×7 median in the run-averaged image. ` +
    `Δg2 is the exact drop in that pixel's ROI g2 when it alone is removed.`;
  const rb = document.getElementById('roibody');
  for (const r of (M.rois || [])){
    const tr = document.createElement('tr');
    const bad = r.nhot > 0;
    tr.innerHTML = `<td>${r.dq}</td><td>${r.q.toFixed(5)}</td>` +
      `<td>${r.nhot || ''}</td><td${bad ? ' style="color:#ff5252"' : ''}>${r.g2.toFixed(3)}</td>` +
      `<td>${r.g2c.toFixed(3)}</td>`;
    tr.addEventListener('click', () => {
      const k = HOT.findIndex(h => h.dq === r.dq); if (k >= 0) gotoHot(k);
    });
    rb.appendChild(tr);
  }
}
function gotoHot(i){
  selHot = i;
  const h = HOT[i], r = stage.getBoundingClientRect();
  z = Math.max(z, 8);
  const yy = flipY ? (NY - 1 - h.y) : h.y;
  tx = r.width / 2 - (h.x + 0.5) * z; ty = r.height / 2 - (yy + 0.5) * z;
  for (const [k, tr] of [...document.getElementById('hotbody').children].entries())
    tr.classList.toggle('sel', k === i);
  blit();
}
document.getElementById('tflag').addEventListener('click', e => {
  showFlags = !showFlags; e.target.textContent = 'markers: ' + (showFlags ? 'on' : 'off'); drawFlags();
});
document.getElementById('tlabel').addEventListener('click', e => {
  showLabels = !showLabels; e.target.textContent = 'labels: ' + (showLabels ? 'on' : 'off'); drawFlags();
});
document.getElementById('prev').addEventListener('click', () =>
  HOT.length && gotoHot((selHot - 1 + HOT.length) % HOT.length));
document.getElementById('next').addEventListener('click', () =>
  HOT.length && gotoHot((selHot + 1) % HOT.length));

cmapSel.addEventListener('change', () => { drawCbar(); renderImage(); });
document.getElementById('scale').addEventListener('change', () => { drawCbar(); renderImage(); });
document.getElementById('rev').addEventListener('click', () => { reversed = !reversed; drawCbar(); renderImage(); });
document.getElementById('flipy').addEventListener('click', () => { flipY = !flipY; renderImage(); });
document.getElementById('fit').addEventListener('click', fitView);
document.getElementById('png').addEventListener('click', () => {
  const a = document.createElement('a');
  a.download = `unnormG2_t${M.i0}_q${M.i1}.png`;
  a.href = off.toDataURL('image/png'); a.click();
});
addEventListener('resize', () => { resize(); });

buildHot();
resize(); fitView(); refresh(false);
</script></body></html>
"""


if __name__ == "__main__":
    main()
