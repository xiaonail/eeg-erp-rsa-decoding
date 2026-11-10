"""step7_specparam.py - Periodic and aperiodic components of the EEG spectrum, face vs car.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import mne
import numpy as np
from fooof import FOOOFGroup
from scipy import io as sio
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROC_DIR = PROJECT_ROOT / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"

FACE_CODES = set(range(1, 41))    # confirm against step1_discover.m output
CAR_CODES = set(range(41, 81))
CONDITIONS = ("face", "car")

TMIN, TMAX = 0.0, 0.8             # post-stimulus window (s)
FREQ_RANGE = (3, 40)              # fitting range (Hz)
BANDWIDTH = 2.5                   # multitaper bandwidth (Hz), the minimum for a 0.8 s window
CHANNELS = ("PO8", "P8", "PO7", "P7")   # preferred channel and fallbacks, as in step 4
MIN_TRIALS = 20


def stimulus_codes(epochs):
    """Returns original ERP CORE stimulus code per epoch"""
    code_of = {internal: int(label.split("/")[0])
               for label, internal in epochs.event_id.items()}
    return np.array([code_of[i] for i in epochs.events[:, 2]])


def condition_spectra(set_path):
    """Returns frequencies, channel names and mean power spectra per condition"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        epochs = mne.io.read_epochs_eeglab(str(set_path), verbose="ERROR")
    epochs.pick("eeg")

    codes = stimulus_codes(epochs)
    spectra, n_trials = {}, {}
    for condition, code_set in zip(CONDITIONS, (FACE_CODES, CAR_CODES)):
        selected = epochs[np.isin(codes, list(code_set))]
        n_trials[condition] = len(selected)
        psd = selected.compute_psd(method="multitaper", fmin=FREQ_RANGE[0], fmax=FREQ_RANGE[1],
                                   tmin=TMIN, tmax=TMAX, bandwidth=BANDWIDTH, verbose="ERROR")
        power, freqs = psd.get_data(return_freqs=True)
        spectra[condition] = power.mean(axis=0)
    return freqs, epochs.ch_names, spectra, n_trials


def fit(freqs, spectra):
    """Returns aperiodic offset, exponent and fit R^2 per channel"""
    model = FOOOFGroup(peak_width_limits=[2, 8], max_n_peaks=4, aperiodic_mode="fixed",
                       verbose=False)
    model.fit(freqs, spectra, list(FREQ_RANGE))
    aperiodic = model.get_params("aperiodic_params")
    return aperiodic[:, 0], aperiodic[:, 1], model.get_params("r_squared")


def main():
    """Fits all subjects and writes figures and results"""
    set_files = sorted(PROC_DIR.glob("*_clean.set"))
    if not set_files:
        raise FileNotFoundError(f"No *_clean.set files in {PROC_DIR}. Run steps 0-2 first.")

    subjects, offset, exponent, r_squared, log_power = [], [], [], [], []
    for i, path in enumerate(set_files, 1):
        subject = path.stem.replace("_clean", "")
        freqs, ch_names, spectra, n_trials = condition_spectra(path)
        if min(n_trials.values()) < MIN_TRIALS:
            print(f"[{i}/{len(set_files)}] {subject}: skipped ({n_trials})")
            continue

        fits = [fit(freqs, spectra[c]) for c in CONDITIONS]
        subjects.append(subject)
        offset.append([f[0] for f in fits])
        exponent.append([f[1] for f in fits])
        r_squared.append([f[2] for f in fits])
        log_power.append([np.log10(spectra[c]) for c in CONDITIONS])
        print(f"[{i}/{len(set_files)}] {subject}: mean exponent "
              f"face={fits[0][1].mean():.2f}, car={fits[1][1].mean():.2f}")

    # Arrays are subject x condition x channel (x frequency)
    offset, exponent = np.array(offset), np.array(exponent)
    r_squared, log_power = np.array(r_squared), np.array(log_power)
    n = len(subjects)

    channel = next(c for c in CHANNELS if c in ch_names)
    ch = ch_names.index(channel)
    t_stat, p_val = stats.ttest_rel(exponent[:, 0, ch], exponent[:, 1, ch])
    print(f"\nN={n}  exponent at {channel}: face={exponent[:, 0, ch].mean():.2f}, "
          f"car={exponent[:, 1, ch].mean():.2f}, t({n - 1})={t_stat:+.2f}, p={p_val:.4f}")
    print(f"Median fit R^2 = {np.median(r_squared):.3f}")

    # Figure 10: spectra with aperiodic fits, and paired exponents
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    for c, (condition, colour) in enumerate(zip(CONDITIONS, ("k", "0.55"))):
        ap_fit = offset[:, c, ch].mean() - exponent[:, c, ch].mean() * np.log10(freqs)
        ax1.plot(freqs, log_power[:, c, ch].mean(0), color=colour, lw=2, label=condition)
        ax1.plot(freqs, ap_fit, color=colour, ls="--", lw=1)
    ax1.set(xscale="log", xticks=[4, 10, 20, 40], xticklabels=["4", "10", "20", "40"],
            xlabel="Frequency (Hz)", ylabel="log10 power",
            title=f"Post-stimulus spectra at {channel}, dashed = aperiodic fit")
    ax1.xaxis.set_minor_formatter(NullFormatter())
    ax1.legend(frameon=False)

    ax2.plot([0, 1], exponent[:, :, ch].T, color="0.6", marker="o")
    ax2.plot([0, 1], exponent[:, :, ch].mean(0), color="k", marker="o", lw=2.5)
    ax2.set(xticks=[0, 1], xticklabels=CONDITIONS, xlim=(-0.3, 1.3), ylabel="Aperiodic exponent",
            title=f"Face vs car, t({n - 1})={t_stat:.2f}, " + (f"p={p_val:.3f}" if p_val >= 1e-3 else "p<0.001"))
    fig.suptitle(f"Spectral parameterisation (N={n})")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig10_specparam.png", dpi=150)
    plt.close(fig)

    sio.savemat(RESULTS_DIR / "specparam_results.mat", {
        "subjects": np.array(subjects, dtype=object), "conditions": np.array(CONDITIONS, dtype=object),
        "ch_names": np.array(ch_names, dtype=object), "freqs": freqs,
        "offset": offset, "exponent": exponent, "r_squared": r_squared,
        "channel": channel, "t_stat": t_stat, "p_val": p_val,
        "window_s": np.array([TMIN, TMAX]), "freq_range": np.array(FREQ_RANGE),
    })

    print("\nSpectral parameterisation complete. Next: code/step8_summary.m")


if __name__ == "__main__":
    main()
