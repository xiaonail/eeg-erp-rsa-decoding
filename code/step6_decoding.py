"""step6_decoding.py - Multivariate decoding of face vs car from EEG patterns.

Trains a regularised logistic classifier on each time point's channel pattern
and scores it by cross-validated AUC. A second pass trains at every time point
and tests at every other, producing a temporal generalisation matrix that shows
whether the face/car distinction is coded briefly or held across time.

Where step5_rsa.m tests whether condition-level geometry matches a category
model, this step asks whether individual trials can be told apart.
"""

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import numpy as np
from mne.decoding import GeneralizingEstimator, SlidingEstimator, cross_val_multiscore
from mne.stats import permutation_cluster_1samp_test
from scipy import io as sio
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROC_DIR = PROJECT_ROOT / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"

FACE_CODES = set(range(1, 41))    # confirm against step1_discover.m output
CAR_CODES = set(range(41, 81))

N_FOLDS = 5
RANDOM_STATE = 42
C_PARAM = 0.1                     # inverse regularisation strength
MIN_TRIALS = 20                   # per class, after artefact rejection
SFREQ_DECODE = 512                # generalisation fits n_times^2 models, so it
SFREQ_GEN = 128                   # runs coarser than the sliding analysis
CHANCE = 0.5
ALPHA = 0.05
N_PERMUTATIONS = 10000


def cluster_test(auc, alpha=ALPHA):
    """Returns significance mask over time and cluster p-values"""
    _, clusters, cluster_pv, _ = permutation_cluster_1samp_test(
        auc - CHANCE, n_permutations=N_PERMUTATIONS, tail=1,
        seed=RANDOM_STATE, out_type="mask", verbose=False)

    sig = np.zeros(auc.shape[1], dtype=bool)
    for cluster, pv in zip(clusters, cluster_pv):
        if pv < alpha:
            sig[cluster] = True   # clusters index directly as tuples
    return sig, np.sort(cluster_pv)


def stimulus_codes(epochs):
    """Returns original ERP CORE stimulus code per epoch"""
    code_of = {internal: int(label.split("/")[0])
               for label, internal in epochs.event_id.items()}
    return np.array([code_of[i] for i in epochs.events[:, 2]])


def load_subject(set_path):
    """Returns face/car epochs and labels, 0 = face and 1 = car"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        epochs = mne.io.read_epochs_eeglab(str(set_path), verbose="ERROR")

    epochs.pick("eeg")  # a stray EOG channel would give the classifier a free ride

    codes = stimulus_codes(epochs)
    is_car = np.isin(codes, list(CAR_CODES))
    keep = np.isin(codes, list(FACE_CODES)) | is_car
    if not keep.any():
        raise ValueError(
            f"No face/car epochs in {set_path.name}; found {sorted(set(codes))[:8]}"
        )

    return epochs[keep], is_car[keep].astype(int)


def classifier():
    """Returns standardised logistic regression refit inside every fold"""
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(C=C_PARAM, solver="liblinear", max_iter=2000,
                           class_weight="balanced", random_state=RANDOM_STATE),
    )


def score(estimator, epochs, y, sfreq, cv):
    """Returns cross-validated AUC averaged over folds"""
    X = epochs.copy().resample(sfreq).get_data(copy=True)
    est = estimator(classifier(), scoring="roc_auc", verbose="ERROR")
    return cross_val_multiscore(est, X, y, cv=cv, verbose="ERROR").mean(axis=0)


def main():
    """Runs decoding for all subjects and writes figures and results"""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    set_files = sorted(PROC_DIR.glob("*_clean.set"))
    if not set_files:
        raise FileNotFoundError(
            f"No *_clean.set files in {PROC_DIR}\n"
            "Run steps 0-2 in MATLAB first; step2_preprocess.m writes them there."
        )
    print(f"Reading {len(set_files)} cleaned files from {PROC_DIR}\n")

    cv = StratifiedKFold(N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    subjects, auc_time, auc_gen, times, times_gen = [], [], [], None, None

    for i, path in enumerate(set_files, 1):
        subject = path.stem.replace("_clean", "")
        epochs, y = load_subject(path)
        n_face, n_car = (y == 0).sum(), y.sum()

        if min(n_face, n_car) < MIN_TRIALS:
            print(f"[{i}/{len(set_files)}] {subject}: skipped "
                  f"({n_face} face / {n_car} car)")
            continue

        scores = score(SlidingEstimator, epochs, y, SFREQ_DECODE, cv)
        if times is None:
            times = epochs.copy().resample(SFREQ_DECODE).times * 1000
            times_gen = epochs.copy().resample(SFREQ_GEN).times * 1000

        subjects.append(subject)
        auc_time.append(scores)
        auc_gen.append(score(GeneralizingEstimator, epochs, y, SFREQ_GEN, cv))

        peak = np.argmax(scores)
        print(f"[{i}/{len(set_files)}] {subject}: peak AUC={scores[peak]:.3f} "
              f"at {times[peak]:.0f} ms ({n_face} face / {n_car} car)")

    auc_time, auc_gen = np.array(auc_time), np.array(auc_gen)
    n = len(subjects)
    if n < 2:
        raise RuntimeError(f"Only {n} subject(s) passed the trial-count check.")

    # Group inference: is decoding above chance, correcting over time?
    sig, cluster_pv = cluster_test(auc_time)
    grand, sem = auc_time.mean(0), auc_time.std(0, ddof=1) / np.sqrt(n)

    peak = np.argmax(grand)
    post = np.where((times >= 0) & sig)[0]
    onset = times[post[0]] if post.size else np.nan

    n170 = auc_time[:, (times >= 140) & (times <= 200)].mean(1)
    t_n170, p_n170 = stats.ttest_1samp(n170, CHANCE)

    n_sig_clusters = int((cluster_pv < ALPHA).sum())
    print(f"\nN={n}  {n_sig_clusters} significant cluster(s), "
          f"{sig.sum()}/{len(times)} time points (p<{ALPHA}, "
          f"{N_PERMUTATIONS} permutations)")
    if cluster_pv.size:
        print(f"Smallest cluster p = {cluster_pv[0]:.4f}")
    print(f"Peak AUC {grand[peak]:.3f} at {times[peak]:.0f} ms, onset {onset:.0f} ms")
    print(f"N170 window: AUC={n170.mean():.3f}, t({n-1})={t_n170:+.2f}, p={p_n170:.4f}")

    # Figure 8: time-resolved decoding
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.axvspan(140, 200, color="0.5", alpha=0.08, label="N170 window")
    ax.fill_between(times, grand - sem, grand + sem, color="0.85", label="SEM")
    ax.plot(times, grand, "k-", lw=2, label="Mean AUC")
    ax.axhline(CHANCE, color="k", ls="--", lw=1)
    ax.axvline(0, color="k", ls="--", lw=1)

    # Fix limits before drawing markers so they stay inside the axes
    lo, hi = (grand - sem).min(), (grand + sem).max()
    pad = 0.1 * max(hi - lo, 1e-6)
    ax.set_ylim(lo - pad, hi + pad)
    if sig.any():
        ax.plot(times[sig], np.full(sig.sum(), hi + 0.5 * pad), "k.", ms=4,
                label=f"cluster p<{ALPHA}")

    ax.set(xlabel="Time (ms)", ylabel="Cross-validated AUC (+/- SEM)",
           title=f"Face vs car decoding (N={n}, {N_FOLDS}-fold CV)")
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig8_time_resolved_decoding.png", dpi=150)
    plt.close(fig)

    # Figure 9: temporal generalisation, diverging colours centred on chance
    grand_gen = auc_gen.mean(0)
    vmax = np.abs(grand_gen - CHANCE).max()
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.imshow(grand_gen, origin="lower", aspect="auto", cmap="RdBu_r",
                   extent=[times_gen[0], times_gen[-1], times_gen[0], times_gen[-1]],
                   vmin=CHANCE - vmax, vmax=CHANCE + vmax)
    ax.axvline(0, color="k", ls="--", lw=0.8)
    ax.axhline(0, color="k", ls="--", lw=0.8)
    ax.set(xlabel="Testing time (ms)", ylabel="Training time (ms)",
           title=f"Temporal generalisation (N={n})")
    fig.colorbar(im, ax=ax, label="AUC")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "fig9_temporal_generalisation.png", dpi=150)
    plt.close(fig)

    np.savez(RESULTS_DIR / "decoding_results.npz",
             subjects=np.array(subjects), times=times, times_gen=times_gen,
             auc_time=auc_time, auc_gen=auc_gen, grand=grand, sem=sem,
             sig=sig, cluster_pv=cluster_pv, n170=n170, onset=onset,
             n_folds=N_FOLDS, c_param=C_PARAM, n_permutations=N_PERMUTATIONS,
             sfreq_decode=SFREQ_DECODE, sfreq_gen=SFREQ_GEN)

    # step7_summary.m cannot read .npz
    sio.savemat(RESULTS_DIR / "decoding_results.mat", {
        "subjects": np.array(subjects, dtype=object),
        "times": times, "grand": grand, "sem": sem,
        "sig": sig, "cluster_pv": cluster_pv,
        "n170": n170, "onset": onset,
        "peak_auc": grand[peak], "peak_ms": times[peak],
        "n_folds": N_FOLDS, "c_param": C_PARAM, "sfreq_decode": SFREQ_DECODE,
    })

    print(f"\nDecoding complete. Next: code/step7_summary.m")


if __name__ == "__main__":
    main()
