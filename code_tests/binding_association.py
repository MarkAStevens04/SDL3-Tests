"""Plot binding affinity (Kd / Ki) distributions for compounds across a list
of targets, pulled from ChEMBL.

Note: This script originally tried to plot k_on vs k_off, but ChEMBL's
kon/k_off rows for these targets have NULL standard_value — only equilibrium
constants (Kd, Ki, IC50) are populated. So we plot affinity instead.

Requirements:
    pip install chembl_webresource_client pandas matplotlib
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client

TARGETS: list[tuple[str, str]] = [
    ("CHEMBL206",  "ESR1 (estrogen receptor α)"),
    ("CHEMBL279",  "VEGFR2"),
    ("CHEMBL301",  "COX-2"),
    ("CHEMBL1862", "ABL1"),
    ("CHEMBL1827", "PI3K-α"),
    ("CHEMBL2842", "mTOR"),
]

AFFINITY_TYPES = ["Kd", "Ki"]

DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(exist_ok=True)
AFFINITY_CSV = DATA_DIR / "chembl_affinity_raw.csv"
CLEAN_CSV = DATA_DIR / "chembl_affinity_clean.csv"
PLOT_PATH = DATA_DIR / "affinity_distribution.png"

# Conversion factors from each unit string into nanomolar (nM).
UNIT_TO_NM = {
    "nM": 1.0,
    "uM": 1_000.0, "µM": 1_000.0, "microM": 1_000.0,
    "mM": 1_000_000.0,
    "M":  1_000_000_000.0,
    "pM": 1e-3,
    "fM": 1e-6,
}


def fetch_affinity(targets: list[tuple[str, str]]) -> pd.DataFrame:
    """Fetch Kd and Ki activity records for each target from ChEMBL."""
    if AFFINITY_CSV.exists():
        print(f"Loading cached affinity data from {AFFINITY_CSV}")
        return pd.read_csv(AFFINITY_CSV)

    activity = new_client.activity
    frames = []
    for chembl_id, label in targets:
        print(f"Querying ChEMBL affinity (Kd/Ki) for {label} ({chembl_id})...")
        records = activity.filter(
            target_chembl_id=chembl_id,
            standard_type__in=AFFINITY_TYPES,
            standard_value__isnull=False,
        ).only(
            "molecule_chembl_id",
            "target_chembl_id",
            "standard_type",
            "standard_value",
            "standard_units",
            "standard_relation",
        )
        df = pd.DataFrame.from_records(records)
        if df.empty:
            print(f"  no records returned for {chembl_id}")
            continue
        df["target_label"] = label
        frames.append(df)
        print(f"  retrieved {len(df)} records")

    if not frames:
        raise RuntimeError("ChEMBL returned no affinity records for any target.")

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(AFFINITY_CSV, index=False)
    return combined


def clean_affinity(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize to nM and compute pAffinity = -log10(M)."""
    df = df.copy()
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value"])
    df = df[df["standard_value"] > 0]
    df = df[df["standard_units"].isin(UNIT_TO_NM)]

    # Keep only exact relations ('=' / '~'); drop censored values like '>'/'<'.
    df = df[df["standard_relation"].isin(["=", "~", None]) | df["standard_relation"].isna()]

    df["value_nM"] = df["standard_value"] * df["standard_units"].map(UNIT_TO_NM)
    df["pAffinity"] = -np.log10(df["value_nM"] * 1e-9)

    # Per (target, molecule, type), collapse replicates with the median.
    clean = (
        df.groupby(
            ["target_chembl_id", "target_label", "molecule_chembl_id", "standard_type"],
            as_index=False,
        )["pAffinity"].median()
    )
    clean.to_csv(CLEAN_CSV, index=False)
    print(f"Cleaned {len(clean)} (target, molecule, type) affinity measurements.")
    return clean


def plot_affinity(clean: pd.DataFrame, out: Path = PLOT_PATH) -> Path:
    """Box + strip plot of pAffinity per target, split by Kd vs Ki."""
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(10.0, 6.5))

    targets = sorted(clean["target_label"].unique())
    types = ["Kd", "Ki"]
    cmap = {"Kd": "#1f77b4", "Ki": "#d62728"}
    width = 0.36

    positions_by_type = {"Kd": [], "Ki": []}
    for i, label in enumerate(targets):
        for j, t in enumerate(types):
            sub = clean[(clean["target_label"] == label) &
                        (clean["standard_type"] == t)]["pAffinity"].values
            if len(sub) == 0:
                continue
            pos = i + (j - 0.5) * width
            positions_by_type[t].append(pos)
            bp = ax.boxplot(
                sub, positions=[pos], widths=width * 0.85,
                patch_artist=True, showfliers=False,
            )
            for box in bp["boxes"]:
                box.set(facecolor=cmap[t], alpha=0.35, edgecolor=cmap[t])
            for med in bp["medians"]:
                med.set(color=cmap[t], linewidth=1.6)
            jitter = (np.random.RandomState(i * 10 + j).rand(len(sub)) - 0.5) * width * 0.6
            ax.scatter(
                np.full_like(sub, pos, dtype=float) + jitter, sub,
                s=12, alpha=0.55, color=cmap[t], edgecolor="white", linewidth=0.3,
            )

    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(targets, rotation=20, ha="right")
    ax.set_ylabel(r"pAffinity = $-\log_{10}(K\ [\mathrm{M}])$", fontsize=12)
    ax.set_title("Binding affinity distributions across targets (ChEMBL)", fontsize=13)

    counts = (clean.groupby(["target_label", "standard_type"])
                   .size().unstack(fill_value=0)
                   .reindex(targets))
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=cmap[t], alpha=0.4, edgecolor=cmap[t],
                      label=f"{t}  (n={counts[t].sum() if t in counts else 0})")
        for t in types
    ]
    ax.legend(handles=legend_handles, loc="best", frameon=True, title="Measurement")

    for ref, tag in [(9, "1 nM"), (6, "1 µM"), (3, "1 mM")]:
        ax.axhline(ref, color="0.7", linestyle="--", linewidth=0.7, zorder=0)
        ax.text(len(targets) - 0.5, ref, f" K = {tag}",
                color="0.45", fontsize=8, va="bottom", ha="right")

    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved plot to {out}")
    return out


def main() -> None:
    raw = fetch_affinity(TARGETS)
    clean = clean_affinity(raw)
    if clean.empty:
        print("No affinity measurements survived cleaning — nothing to plot.")
        return
    plot_affinity(clean)


if __name__ == "__main__":
    main()
