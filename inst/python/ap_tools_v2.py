import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# reuse your stim_window_from_command from wherever you keep it
# from .whatever import stim_window_from_command

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

def detect_spikes_dvdt(
    *,
    t_s: np.ndarray,
    v_mV: np.ndarray,
    stim_start_s: float | None = None,
    stim_end_s: float | None = None,
    only_stim: bool = True,
    smooth_ms: float = 0.2,
    dvdt_thr: float = 25.0,          # mV/ms (used to estimate AP_begin; can be None to disable)
    refractory_ms: float = 2.0,
    peak_search_ms: float = 2.0,
    return_ap_begin_times: bool = False,
    # --- new / tunable, no hardcoded cutoffs ---
    min_height_mV: float | None = None,       # if None -> derived from noise (RMS)
    prominence_mV: float | None = None,       # if None -> derived from noise (RMS)
    noise_win_s: float = 0.05,                # window used to estimate noise (prefer baseline just before stim)
    noise_mult_height: float = 6.0,           # height >= baseline + mult*RMS
    noise_mult_prom: float = 5.0,             # prominence >= mult*RMS
    min_height_floor_mV: float = 6.0,         # don't let derived thresholds go below this
    prominence_floor_mV: float = 4.0,         # don't let derived thresholds go below this
    min_peak_voltage_mV: float | None = -20.0,
):
    """
    Spike detector using peak prominence + (optional) dv/dt for AP-begin estimation.

    Returns:
      peak_times_s, peak_indices_full
      optionally ap_begin_times_s (estimated as first dv/dt crossing before peak)

    Notes:
      - Uses find_peaks on the *analysis segment* (stim-only or full), then maps indices back to FULL trace.
      - Height/prominence are either user-specified or derived from baseline noise (RMS).
      - No hardcoded 'height=8' / 'prominence=6' etc.
    """
    import numpy as np
    from scipy.signal import find_peaks

    t_s = np.asarray(t_s, float)
    v_mV = np.asarray(v_mV, float)

    if t_s.size < 5 or v_mV.size != t_s.size:
        if return_ap_begin_times:
            return np.array([]), np.array([], dtype=int), np.array([])
        return np.array([]), np.array([], dtype=int)

    # -----------------------------
    # choose analysis segment (stim-only) BUT keep mapping to full indices
    # -----------------------------
    t_full = t_s
    v_full = v_mV

    if (
        only_stim
        and stim_start_s is not None
        and stim_end_s is not None
        and np.isfinite(stim_start_s)
        and np.isfinite(stim_end_s)
        and stim_end_s > stim_start_s
    ):
        m = (t_full >= stim_start_s) & (t_full <= stim_end_s)
        if np.count_nonzero(m) < 5:
            if return_ap_begin_times:
                return np.array([]), np.array([], dtype=int), np.array([])
            return np.array([]), np.array([], dtype=int)
    else:
        m = np.ones_like(t_full, dtype=bool)

    idx_full = np.flatnonzero(m)   # seg-index -> full-index mapping
    t = t_full[m]
    v = v_full[m]

    # dt on analysis segment
    dt = float(np.median(np.diff(t)))
    if not np.isfinite(dt) or dt <= 0:
        if return_ap_begin_times:
            return np.array([]), np.array([], dtype=int), np.array([])
        return np.array([]), np.array([], dtype=int)

    # -----------------------------
    # smoothing for dv/dt + (optionally) for peak finding
    # -----------------------------
    if smooth_ms and smooth_ms > 0:
        win = int(round((smooth_ms / 1000.0) / dt))
        win = max(1, win)
        if win > 1:
            k = np.ones(win, dtype=float) / float(win)
            v_s = np.convolve(v, k, mode="same")
        else:
            v_s = v
    else:
        v_s = v

    # dv/dt in mV/ms (computed on smoothed trace)
    dvdt = np.gradient(v_s, t) / 1000.0

    # -----------------------------
    # baseline + noise estimate (RMS)
    # -----------------------------
    def _estimate_noise_rms() -> tuple[float, float]:
        """
        Returns (baseline_mV, rms_mV) using a window:
          - if stim times exist: the last noise_win_s BEFORE stim_start
          - else: first noise_win_s of the analysis segment
        """
        if (
            stim_start_s is not None
            and np.isfinite(stim_start_s)
            and np.any(t_full < stim_start_s)
        ):
            t1 = float(stim_start_s)
            t0 = max(float(t_full[0]), t1 - float(noise_win_s))
            mm = (t_full >= t0) & (t_full <= t1)
            vv = v_full[mm]
            if vv.size < 10:
                # fallback to segment start
                t0 = float(t[0])
                t1 = float(min(t[-1], t[0] + float(noise_win_s)))
                mm2 = (t >= t0) & (t <= t1)
                vv = v_s[mm2]
        else:
            t0 = float(t[0])
            t1 = float(min(t[-1], t[0] + float(noise_win_s)))
            mm2 = (t >= t0) & (t <= t1)
            vv = v_s[mm2]

        if vv.size < 5:
            baseline = float(np.median(v_s[: max(5, int(0.1 * v_s.size))]))
            rms = float(np.sqrt(np.mean((v_s - baseline) ** 2)))
            return baseline, rms

        baseline = float(np.median(vv))
        rms = float(np.sqrt(np.mean((vv - baseline) ** 2)))
        return baseline, rms

    baseline_mV, noise_rms_mV = _estimate_noise_rms()
    if not np.isfinite(noise_rms_mV) or noise_rms_mV <= 0:
        # ultra-safe fallback
        noise_rms_mV = float(np.nanstd(v_s)) if np.isfinite(np.nanstd(v_s)) and np.nanstd(v_s) > 0 else 1.0

    # derived thresholds if not provided
    if min_height_mV is None:
        min_height_mV = max(float(min_height_floor_mV), float(noise_mult_height) * float(noise_rms_mV))
    if prominence_mV is None:
        prominence_mV = max(float(prominence_floor_mV), float(noise_mult_prom) * float(noise_rms_mV))

    # -----------------------------
    # peak finding on the ANALYSIS segment (indices are segment-local)
    # -----------------------------
    # refractory converted to min distance in samples
    refr_s = float(refractory_ms) / 1000.0
    min_dist = max(1, int(round(refr_s / dt)))

    peaks_seg, props = find_peaks(
        v_s,
        height=baseline_mV + float(min_height_mV),
        prominence=float(prominence_mV),
        distance=min_dist,
    )

    if peaks_seg.size == 0:
        if return_ap_begin_times:
            return np.array([]), np.array([], dtype=int), np.array([])
        return np.array([]), np.array([], dtype=int)

    # -------------------------------------------------
    # ADD THIS BLOCK *RIGHT HERE* (filter by peak V)
    # -------------------------------------------------
    if min_peak_voltage_mV is not None:
        keep = v_s[peaks_seg] >= float(min_peak_voltage_mV)
        peaks_seg = peaks_seg[keep]
        if peaks_seg.size == 0:
            if return_ap_begin_times:
                return np.array([]), np.array([], dtype=int), np.array([])
            return np.array([]), np.array([], dtype=int)
    # -------------------------------------------------

    # map segment peaks -> full indices
    peaks_full = idx_full[peaks_seg]

    # -----------------------------
    # AP begin estimation (optional): first dv/dt crossing before peak
    # -----------------------------
    if return_ap_begin_times:
        search_s = float(peak_search_ms) / 1000.0
        search_n = max(3, int(round(search_s / dt)))

        begins_full = np.empty_like(peaks_full)
        for k, p_seg in enumerate(peaks_seg):
            j1 = int(p_seg)
            j0 = max(0, j1 - search_n)

            if dvdt_thr is None or not np.isfinite(dvdt_thr):
                # fallback: just use peak index if dv/dt threshold disabled
                begin_seg = j1
            else:
                # find first upward crossing of dvdt_thr in [j0, j1]
                seg = dvdt[j0 : j1 + 1]
                cross = np.flatnonzero(seg >= float(dvdt_thr))
                begin_seg = (j0 + int(cross[0])) if cross.size else j0

            begins_full[k] = idx_full[begin_seg]

        peak_times_s = t_full[peaks_full]
        ap_begin_times_s = t_full[begins_full]
        return peak_times_s, peaks_full.astype(int), ap_begin_times_s

    peak_times_s = t_full[peaks_full]
    return peak_times_s, peaks_full.astype(int)


def phase_plane_metrics_one_ap(
    t_s,
    v_mV,
    ap_begin_s,
    peak_s,
    window_after_peak_ms=3.0,
):
    t_s = np.asarray(t_s, float)
    v_mV = np.asarray(v_mV, float)

    if not (np.isfinite(ap_begin_s) and np.isfinite(peak_s)):
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    t0 = float(ap_begin_s)
    t1 = float(peak_s) + float(window_after_peak_ms) / 1000.0
    if t1 <= t0:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    m = (t_s >= t0) & (t_s <= t1)
    if np.count_nonzero(m) < 5:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    t_seg = t_s[m]
    v_seg = v_mV[m]
    dvdt = np.gradient(v_seg, t_seg) / 1000.0  # mV/ms

    i_max = int(np.nanargmax(dvdt))
    i_min = int(np.nanargmin(dvdt))

    return (
        float(dvdt[i_max]), float(v_seg[i_max]), float(t_seg[i_max]),
        float(dvdt[i_min]), float(v_seg[i_min]), float(t_seg[i_min]),
    )


def efel_spikes_for_well(
    df_well,        # polars DF
    v_is_volts=True,
    cfg=None,       # pass your EfelConfig
):
    """
    COMPATIBILITY NAME KEPT.
    Now returns:
      spikes_df: per spike with peak + begin + dvdt extrema points
      sweep_df: per sweep stim window + spike_count
    """
    pdf = df_well.to_pandas()

    spikes_rows = []
    sweep_rows = []

    # defaults if cfg not provided
    smooth_ms = getattr(cfg, "smooth_ms", 0.2) if cfg is not None else 0.2
    dvdt_thr = getattr(cfg, "dvdt_thr_mV_per_ms", 25.0) if cfg is not None else 25.0
    min_height = getattr(cfg, "min_height_mV", 25.0) if cfg is not None else 25.0
    refractory_ms = getattr(cfg, "refractory_ms", 2.0) if cfg is not None else 2.0
    peak_search_ms = getattr(cfg, "peak_search_ms", 2.0) if cfg is not None else 2.0
    stim_ignore_ms = getattr(cfg, "stim_start_ignore_ms", 2.0) if cfg is not None else 2.0
    phase_after_ms = getattr(cfg, "phase_window_after_peak_ms", 3.0) if cfg is not None else 3.0

    for sweep, g in pdf.groupby("sweep", sort=True):
        t_s = g["t_s"].to_numpy(float)
        stim = g["stimulus"].to_numpy(float)
        rec = g["recorded"].to_numpy(float)
        v_mV = rec * 1000.0 if v_is_volts else rec

        t0, t1, amp = stim_window_from_command(t_s, stim)

        t0_eff = t0 + stim_ignore_ms / 1000.0 if np.isfinite(t0) else np.nan

        peak_times_s, peak_idx, ap_begin_times_s = detect_spikes_dvdt(
            t_s=t_s,
            v_mV=v_mV,
            stim_start_s=t0_eff,
            stim_end_s=t1,
            only_stim=True,
            smooth_ms=smooth_ms,
            dvdt_thr=dvdt_thr,
            min_height_mV=min_height,
            refractory_ms=refractory_ms,
            peak_search_ms=peak_search_ms,
            return_ap_begin_times=True,
        )

        spike_count = int(len(peak_times_s))

        for st, ab in zip(peak_times_s, ap_begin_times_s):
            dvdt_max, dvdt_max_v, dvdt_max_t, dvdt_min, dvdt_min_v, dvdt_min_t = phase_plane_metrics_one_ap(
                t_s=t_s,
                v_mV=v_mV,
                ap_begin_s=float(ab),
                peak_s=float(st),
                window_after_peak_ms=float(phase_after_ms),
            )

            spikes_rows.append({
                "sweep": int(sweep),
                "spike_time_s": float(st),          # peak time
                "ap_begin_time_s": float(ab),
                "dvdt_max_time_s": float(dvdt_max_t),
                "dvdt_max_v_mV": float(dvdt_max_v),
                "dvdt_min_time_s": float(dvdt_min_t),
                "dvdt_min_v_mV": float(dvdt_min_v),
            })

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


# -----------------------------
# Plotting (NO plt.show() here)
# -----------------------------

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
            plt.scatter(sp, ysp, s=12)  # peak

        # mark dv/dt extrema points if present
        if "dvdt_max_time_s" in spikes_df.columns:
            tmax = spikes_df.loc[spikes_df["sweep"] == sweep, "dvdt_max_time_s"].to_numpy(float)
            if tmax.size:
                ymax = np.interp(tmax, t_s, v_mV)
                plt.scatter(tmax, ymax, s=16, marker="^")  # dvdt max

        if "dvdt_min_time_s" in spikes_df.columns:
            tmin = spikes_df.loc[spikes_df["sweep"] == sweep, "dvdt_min_time_s"].to_numpy(float)
            if tmin.size:
                ymin = np.interp(tmin, t_s, v_mV)
                plt.scatter(tmin, ymin, s=16, marker="v")  # dvdt min

    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (mV)")
    plt.title("All sweeps (overlay) with spike peaks + dv/dt extrema")


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
            plt.scatter(sp, ysp, s=14)

        if "dvdt_max_time_s" in spikes_df.columns:
            tmax = spikes_df.loc[spikes_df["sweep"] == sweep, "dvdt_max_time_s"].to_numpy(float)
            if tmax.size:
                ymax = np.interp(tmax, t_s, y)
                plt.scatter(tmax, ymax, s=16, marker="^")

        if "dvdt_min_time_s" in spikes_df.columns:
            tmin = spikes_df.loc[spikes_df["sweep"] == sweep, "dvdt_min_time_s"].to_numpy(float)
            if tmin.size:
                ymin = np.interp(tmin, t_s, y)
                plt.scatter(tmin, ymin, s=16, marker="v")

    plt.xlabel("Time (s)")
    plt.ylabel("Voltage + offset (mV)")
    plt.title("All sweeps (stacked) with spike peaks + dv/dt extrema")


def phase_plane_many_sweeps(
    df_well,
    sweeps=None,
    v_is_volts=True,
    smooth_ms=0.3,
    stim_window_df=None,
    only_stim=True,
    max_points_per_sweep=8000,
    alpha=0.35
):
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

        if only_stim and stim_window_df is not None:
            row = stim_window_df[stim_window_df["sweep"] == s]
            if len(row):
                t0 = float(row["stim_start_s"].iloc[0])
                t1 = float(row["stim_end_s"].iloc[0])
                if np.isfinite(t0) and np.isfinite(t1):
                    m = (t_s >= t0) & (t_s <= t1)
                    t_s = t_s[m]
                    v = v[m]

        v_mV = v * 1000.0 if v_is_volts else v
        if len(v_mV) < 5:
            continue
        dt = np.median(np.diff(t_s))
        if not np.isfinite(dt) or dt <= 0:
            continue

        if smooth_ms and smooth_ms > 0:
            win = int(round((smooth_ms / 1000.0) / dt))
            win = max(1, win)
            if win > 1:
                kernel = np.ones(win) / win
                v_mV = np.convolve(v_mV, kernel, mode="same")

        dvdt = np.gradient(v_mV, t_s) / 1000.0

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