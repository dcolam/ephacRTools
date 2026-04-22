from pathlib import Path
import re
import pandas as pd
import numpy as np

import pyarrow as pa
import pyarrow.parquet as pq


# -------- helpers --------

SWEEP_SUFFIX = re.compile(r"_(\d{3})\.csv$", re.IGNORECASE)

from pathlib import Path
import polars as pl

def load_well_from_parquet(
    parquet_path,
    plate_id,
    well_id,
    cols=("plate_id", "well_id", "sweep", "t_s", "recorded", "stimulus"),
):
    parquet_path = str(parquet_path)

    df = (
        pl.scan_parquet(parquet_path)
        .filter((pl.col("plate_id") == plate_id) & (pl.col("well_id") == well_id))
        .select([c for c in cols if c is not None])
        .collect()
        .sort(["sweep", "t_s"])
    )
    return df

def parse_filename(path: Path):
    """
    Example filename:
      251223_VI 10pA_16.30.19_A01_22W16255.csv
    Returns dict with fields.
    """
    stem = path.stem
    parts = stem.split("_")
    # robust-ish: last two are well and plate
    plate_id = parts[-1]
    well_id = parts[-2]
    protocol = "_".join(parts[1:-2])  # everything between date and well
    date = parts[0]
    return {"date": date, "protocol": protocol, "well_id": well_id, "plate_id": plate_id}


def read_sweep_times_first_line(path: Path):
    """
    First line is like:
      ;Sweep Time (s);;;603.000E-3;;12.566E+0;;...
    We extract the numeric tokens in order.
    """
    first = path.read_text(errors="ignore").splitlines()[0]
    toks = [t.strip() for t in first.split(";") if t.strip()]
    sweep_times = []
    for t in toks:
        try:
            sweep_times.append(float(t))
        except ValueError:
            pass
    return np.array(sweep_times, dtype=float)


def read_well_csv_to_long(path: Path):
    meta = parse_filename(path)
    sweep_times = read_sweep_times_first_line(path)

    # main table starts at line 2 (skip first line)
    df = pd.read_csv(path, sep=";", skiprows=1, engine="python")

    # time axis (seconds)
    if "Sample Time (us)" in df.columns:
        t_s = df["Sample Time (us)"].astype(float).to_numpy() / 1e6
    else:
        # fallback: sample index
        t_s = np.arange(len(df), dtype=float)

    # identify sweep columns
    raw_cols  = [c for c in df.columns if c.startswith("Raw Data")]
    stim_cols = [c for c in df.columns if c.startswith("Stimulus")]

    n_sweeps = min(len(raw_cols), len(stim_cols))
    raw_cols = raw_cols[:n_sweeps]
    stim_cols = stim_cols[:n_sweeps]

    if len(sweep_times) < n_sweeps:
        # if header had fewer numeric tokens, pad
        sweep_times = np.concatenate([sweep_times, np.full(n_sweeps - len(sweep_times), np.nan)])
    sweep_times = sweep_times[:n_sweeps]

    # build long blocks sweep-by-sweep (memory safe-ish per well)
    blocks = []
    for i in range(n_sweeps):
        blocks.append(pd.DataFrame({
            "date": meta["date"],
            "protocol": meta["protocol"],
            "plate_id": meta["plate_id"],
            "well_id": meta["well_id"],
            "sweep": i + 1,
            "sweep_time_s": sweep_times[i],
            "t_s": t_s,
            "recorded": df[raw_cols[i]].astype(float).to_numpy(),
            "stimulus": df[stim_cols[i]].astype(float).to_numpy(),
        }))

    long_df = pd.concat(blocks, ignore_index=True)
    return long_df


# -------- main: convert folder of 384 wells into ONE parquet --------

def folder_csvs_to_single_parquet(input_dir, output_parquet):
    input_dir = Path(input_dir)
    csvs = sorted(input_dir.glob("*.csv"))
    if not csvs:
        raise FileNotFoundError(f"No CSVs found in {input_dir}")

    writer = None
    outpath = Path(output_parquet)

    for i, f in enumerate(csvs, start=1):
        long_df = read_well_csv_to_long(f)
        table = pa.Table.from_pandas(long_df, preserve_index=False)

        if writer is None:
            writer = pq.ParquetWriter(outpath, table.schema, compression="zstd")

        writer.write_table(table)
        if i % 20 == 0:
            print(f"Wrote {i}/{len(csvs)} wells...")

    if writer is not None:
        writer.close()
    print("Done:", outpath)


import numpy as np
import pandas as pd
import efel

def stim_window_from_command(t_s, stim, frac=0.2, min_dur_s=0.02):
    stim = np.asarray(stim, float)
    t_s = np.asarray(t_s, float)

    n0 = max(1, int(0.1 * len(stim)))
    base = np.median(stim[:n0])
    dev = stim - base

    topk = max(10, len(dev) // 50)
    amp = float(np.median(dev[np.argsort(np.abs(dev))[-topk:]]))
    if not np.isfinite(amp) or abs(amp) < 1e-15:
        return np.nan, np.nan, 0.0

    thr = base + frac * amp
    mask = stim >= thr if amp > 0 else stim <= thr
    idx = np.where(mask)[0]
    if idx.size == 0:
        return np.nan, np.nan, amp

    splits = np.where(np.diff(idx) > 1)[0]
    runs = np.split(idx, splits + 1)
    run = max(runs, key=len)

    t0, t1 = float(t_s[run[0]]), float(t_s[run[-1]])
    if (t1 - t0) < min_dur_s:
        return np.nan, np.nan, amp
    return t0, t1, amp


def efel_spikes_for_well(
    df_well,                  # polars DF
    v_is_volts=True,           # convert recorded to mV if in V
):
    """
    Returns:
      spikes_df: one row per spike (sweep, spike_time_s)
      sweep_df:  one row per sweep with stim metadata and spike count
    """
    pdf = df_well.to_pandas()

    spikes_rows = []
    sweep_rows = []

    # ask eFEL for a few candidate spike-time features; we’ll use whichever comes back
    wanted = [
        "Spikecount",
        "peak_time",          # often available (ms)
        "AP_begin_time",      # fallback (ms)
        "AP_begin_indices",   # fallback (indices)
    ]

    for sweep, g in pdf.groupby("sweep", sort=True):
        t_s = g["t_s"].to_numpy(float)
        stim = g["stimulus"].to_numpy(float)
        rec = g["recorded"].to_numpy(float)

        v_mV = rec * 1000.0 if v_is_volts else rec
        t_ms = t_s * 1000.0

        t0, t1, amp = stim_window_from_command(t_s, stim)

        trace = {
            "T": t_ms,
            "V": v_mV,
            "stim_start": [float(t0 * 1000.0)] if np.isfinite(t0) else [float(t_ms[0])],
            "stim_end":   [float(t1 * 1000.0)] if np.isfinite(t1) else [float(t_ms[-1])],
        }

        # eFEL returns a dict of feature -> list/array/None
        out = efel.getFeatureValues([trace], wanted)[0]

        # spike count
        sc = out.get("Spikecount", None)
        spike_count = int(sc[0]) if isinstance(sc, (list, np.ndarray)) and len(sc) else (0 if sc is None else int(sc))

        # spike times (prefer peak_time in ms, else AP_begin_time in ms, else indices)
        spike_times_s = np.array([], dtype=float)

        if out.get("peak_time", None) is not None and len(out["peak_time"]):
            spike_times_s = np.asarray(out["peak_time"], float) / 1000.0

        elif out.get("AP_begin_time", None) is not None and len(out["AP_begin_time"]):
            spike_times_s = np.asarray(out["AP_begin_time"], float) / 1000.0

        elif out.get("AP_begin_indices", None) is not None and len(out["AP_begin_indices"]):
            idx = np.asarray(out["AP_begin_indices"], int)
            idx = np.clip(idx, 0, len(t_s) - 1)
            spike_times_s = t_s[idx]

        # store spikes
        for st in spike_times_s:
            spikes_rows.append({"sweep": int(sweep), "spike_time_s": float(st)})

        sweep_rows.append({
            "sweep": int(sweep),
            "stim_start_s": float(t0) if np.isfinite(t0) else np.nan,
            "stim_end_s": float(t1) if np.isfinite(t1) else np.nan,
            "stim_amp": float(amp),
            "spike_count": int(spike_count),
        })

    spikes_df = pd.DataFrame(spikes_rows)
    sweep_df = pd.DataFrame(sweep_rows).sort_values("sweep")
    return spikes_df, sweep_df

import matplotlib.pyplot as plt
import numpy as np

def plot_all_sweeps_overlay(df_well, spikes_df, v_is_volts=True, alpha=0.25):
    pdf = df_well.to_pandas()

    plt.figure(figsize=(10, 6))

    for sweep, g in pdf.groupby("sweep", sort=True):
        t_s = g["t_s"].to_numpy(float)
        rec = g["recorded"].to_numpy(float)
        v_mV = rec * 1000.0 if v_is_volts else rec

        plt.plot(t_s, v_mV, alpha=alpha, linewidth=1)

        sp = spikes_df.loc[spikes_df["sweep"] == sweep, "spike_time_s"].to_numpy(float)
        if sp.size:
            ysp = np.interp(sp, t_s, v_mV)
            plt.scatter(sp, ysp, s=10)

    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (mV)")
    plt.title("All sweeps (overlay) with eFEL spike markers")
    plt.show()

def plot_all_sweeps_stacked(df_well, spikes_df, v_is_volts=True, offset_mV=40.0):
    pdf = df_well.to_pandas()

    plt.figure(figsize=(10, 8))

    for i, (sweep, g) in enumerate(pdf.groupby("sweep", sort=True), start=0):
        t_s = g["t_s"].to_numpy(float)
        rec = g["recorded"].to_numpy(float)
        v_mV = rec * 1000.0 if v_is_volts else rec
        y = v_mV + i * offset_mV

        plt.plot(t_s, y, linewidth=0.9)

        sp = spikes_df.loc[spikes_df["sweep"] == sweep, "spike_time_s"].to_numpy(float)
        if sp.size:
            ysp = np.interp(sp, t_s, y)
            plt.scatter(sp, ysp, s=12)

    plt.xlabel("Time (s)")
    plt.ylabel("Voltage + offset (mV)")
    plt.title("All sweeps (stacked) with eFEL spike markers")
    plt.show()
import numpy as np
import matplotlib.pyplot as plt

def phase_plane_one_sweep(
    t_s,
    v,
    v_is_volts=True,
    smooth_ms=0.3,
    max_points=20000,
    title=None,
    spike_times_s=None,   # optional: mark eFEL spike times
):
    """
    Plot dV/dt (mV/ms) vs V (mV) for a single sweep.

    Parameters
    ----------
    t_s : array-like
        Time in seconds.
    v : array-like
        Voltage in Volts or mV.
    v_is_volts : bool
        If True, convert V->mV.
    smooth_ms : float
        Moving-average smoothing window in milliseconds (0 disables).
        Helps reduce derivative noise.
    max_points : int
        Downsample for speed/plot size (None disables).
    spike_times_s : array-like
        Spike times in seconds to highlight on the phase plane.
    """
    t_s = np.asarray(t_s, float)
    v = np.asarray(v, float)

    v_mV = v * 1000.0 if v_is_volts else v

    # estimate dt
    dt = np.median(np.diff(t_s))
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("Bad time axis; cannot compute derivative.")

    # optional smoothing (moving average)
    if smooth_ms and smooth_ms > 0:
        win = int(round((smooth_ms / 1000.0) / dt))
        win = max(1, win)
        if win > 1:
            kernel = np.ones(win) / win
            v_mV_s = np.convolve(v_mV, kernel, mode="same")
        else:
            v_mV_s = v_mV
    else:
        v_mV_s = v_mV

    # derivative: mV/ms
    dvdt = np.gradient(v_mV_s, t_s) / 1000.0  # (mV/s) -> (mV/ms)

    # optional downsample
    if max_points is not None and len(v_mV_s) > max_points:
        idx = np.linspace(0, len(v_mV_s) - 1, max_points).astype(int)
        v_plot = v_mV_s[idx]
        dvdt_plot = dvdt[idx]
        t_plot = t_s[idx]
    else:
        v_plot, dvdt_plot, t_plot = v_mV_s, dvdt, t_s

    plt.figure(figsize=(7, 6))
    plt.plot(v_plot, dvdt_plot, linewidth=1)
    plt.xlabel("V (mV)")
    plt.ylabel("dV/dt (mV/ms)")
    plt.title(title or "Phase-plane: dV/dt vs V")
    plt.axhline(0, linestyle="dashed", linewidth=0.8)
    plt.axvline(0, linestyle="dashed", linewidth=0.8)

    # mark spikes (project to nearest time index)
    if spike_times_s is not None:
        spike_times_s = np.asarray(spike_times_s, float)
        if spike_times_s.size:
            j = np.searchsorted(t_plot, spike_times_s)
            j = np.clip(j, 0, len(t_plot) - 1)
            plt.scatter(v_plot[j], dvdt_plot[j], s=35)

    plt.tight_layout()
    plt.show()
def phase_plane_many_sweeps(
    df_well,                # polars df OR pandas df
    sweeps=None,
    v_is_volts=True,
    smooth_ms=0.3,
    stim_window_df=None,    # optional: sweep_df with stim_start_s/stim_end_s
    only_stim=True,
    max_points_per_sweep=8000,
    alpha=0.35
):
    """
    Overlay phase-plane plots for multiple sweeps from a well.
    """
    # accept polars or pandas
    pdf = df_well.to_pandas() if hasattr(df_well, "to_pandas") else df_well

    if sweeps is None:
        sweeps = sorted(pdf["sweep"].unique())

    plt.figure(figsize=(7, 6))

    for s in sweeps:
        g = pdf[pdf["sweep"] == s]
        if g.empty:
            continue
        t_s = g["t_s"].to_numpy(float)
        v = g["recorded"].to_numpy(float)

        # restrict to stim window if provided
        if only_stim and stim_window_df is not None:
            row = stim_window_df[stim_window_df["sweep"] == s]
            if len(row):
                t0 = float(row["stim_start_s"].iloc[0])
                t1 = float(row["stim_end_s"].iloc[0])
                if np.isfinite(t0) and np.isfinite(t1):
                    m = (t_s >= t0) & (t_s <= t1)
                    t_s = t_s[m]
                    v = v[m]

        # compute phase-plane for this sweep (no individual figure)
        v_mV = v * 1000.0 if v_is_volts else v
        dt = np.median(np.diff(t_s))
        if not np.isfinite(dt) or dt <= 0 or len(v_mV) < 5:
            continue

        if smooth_ms and smooth_ms > 0:
            win = int(round((smooth_ms / 1000.0) / dt))
            win = max(1, win)
            if win > 1:
                kernel = np.ones(win) / win
                v_mV = np.convolve(v_mV, kernel, mode="same")

        dvdt = np.gradient(v_mV, t_s) / 1000.0  # mV/ms

        # downsample per sweep
        if max_points_per_sweep is not None and len(v_mV) > max_points_per_sweep:
            idx = np.linspace(0, len(v_mV) - 1, max_points_per_sweep).astype(int)
            v_mV = v_mV[idx]
            dvdt = dvdt[idx]

        plt.plot(v_mV, dvdt, alpha=alpha, linewidth=1)

    plt.xlabel("V (mV)")
    plt.ylabel("dV/dt (mV/ms)")
    plt.title("Phase-plane overlay: dV/dt vs V")
    plt.axhline(0, linestyle="dashed", linewidth=0.8)
    plt.axvline(0, linestyle="dashed", linewidth=0.8)
    plt.tight_layout()
    plt.show()