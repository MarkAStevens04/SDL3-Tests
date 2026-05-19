"""Download bioactivity data for the estrogen receptor from ChEMBL
and explore the compounds with RDKit.

Requirements:
    pip install chembl_webresource_client rdkit pandas
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Draw, Lipinski

RDLogger.DisableLog("rdApp.*")

# Estrogen receptor alpha (ESR1, human) — ChEMBL target ID
TARGET_CHEMBL_ID = "CHEMBL206"
DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(exist_ok=True)
RAW_CSV = DATA_DIR / "chembl_estrogen_receptor_raw.csv"
CLEAN_CSV = DATA_DIR / "chembl_estrogen_receptor_clean.csv"


def fetch_bioactivities(target_id: str) -> pd.DataFrame:
    """Fetch IC50 bioactivities for the given target from ChEMBL."""
    if RAW_CSV.exists():
        print(f"Loading cached data from {RAW_CSV}")
        return pd.read_csv(RAW_CSV)

    print(f"Querying ChEMBL for bioactivities against {target_id}...")
    activity = new_client.activity
    records = activity.filter(
        target_chembl_id=target_id,
        standard_type="IC50",
    ).only(
        "molecule_chembl_id",
        "canonical_smiles",
        "standard_value",
        "standard_units",
        "standard_relation",
        "assay_type",
        "pchembl_value",
    )
    df = pd.DataFrame.from_records(records)
    print(f"Retrieved {len(df)} activity records.")
    df.to_csv(RAW_CSV, index=False)
    return df


def clean_bioactivities(df: pd.DataFrame) -> pd.DataFrame:
    """Drop rows with missing SMILES or IC50, deduplicate compounds."""
    df = df.dropna(subset=["canonical_smiles", "standard_value"]).copy()
    df = df[df["standard_units"] == "nM"]
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value"])
    df = df[df["standard_value"] > 0]
    # Keep the most potent measurement per molecule
    df = df.sort_values("standard_value").drop_duplicates("molecule_chembl_id")
    df["pIC50"] = -np.log10(df["standard_value"].astype(float) * 1e-9)
    return df.reset_index(drop=True)


def compute_descriptors(smiles: str) -> dict | None:
    """Compute a panel of RDKit descriptors for a single SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "mol_weight": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "h_donors": Lipinski.NumHDonors(mol),
        "h_acceptors": Lipinski.NumHAcceptors(mol),
        "rot_bonds": Lipinski.NumRotatableBonds(mol),
        "tpsa": Descriptors.TPSA(mol),
        "rings": Lipinski.RingCount(mol),
        "aromatic_rings": Lipinski.NumAromaticRings(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
    }


def lipinski_pass(row: pd.Series) -> bool:
    """Lipinski's rule of five: at most one violation."""
    violations = sum([
        row["mol_weight"] > 500,
        row["logp"] > 5,
        row["h_donors"] > 5,
        row["h_acceptors"] > 10,
    ])
    return violations <= 1


def enrich_with_rdkit(df: pd.DataFrame) -> pd.DataFrame:
    print("Computing RDKit descriptors...")
    desc = df["canonical_smiles"].apply(compute_descriptors)
    df = df.join(pd.DataFrame(desc.tolist(), index=df.index))
    df = df.dropna(subset=["mol_weight"])
    df["lipinski_pass"] = df.apply(lipinski_pass, axis=1)
    return df


def draw_top_compounds(df: pd.DataFrame, n: int = 12) -> Path:
    """Render the top-N most potent compounds to a PNG grid."""
    top = df.head(n)
    mols = [Chem.MolFromSmiles(s) for s in top["canonical_smiles"]]
    legends = [
        f"{cid}\nIC50={val:.1f} nM"
        for cid, val in zip(top["molecule_chembl_id"], top["standard_value"])
    ]
    img_path = DATA_DIR / "top_compounds.png"
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(300, 300), legends=legends)
    img.save(img_path)
    return img_path


def summarize(df: pd.DataFrame) -> None:
    print("\n=== Dataset summary ===")
    print(f"Unique compounds:        {len(df)}")
    print(f"Median IC50 (nM):        {df['standard_value'].median():.1f}")
    print(f"Lipinski pass fraction:  {df['lipinski_pass'].mean():.1%}")
    print("\nDescriptor statistics:")
    cols = ["mol_weight", "logp", "h_donors", "h_acceptors", "tpsa", "rot_bonds"]
    print(df[cols].describe().round(2))
    print("\nMost potent compounds:")
    print(
        df.sort_values("standard_value")
          .head(10)[["molecule_chembl_id", "standard_value", "mol_weight", "logp"]]
          .to_string(index=False)
    )


def main() -> None:
    raw = fetch_bioactivities(TARGET_CHEMBL_ID)
    cleaned = clean_bioactivities(raw)
    enriched = enrich_with_rdkit(cleaned)
    enriched.to_csv(CLEAN_CSV, index=False)
    print(f"\nSaved cleaned dataset to {CLEAN_CSV}")
    summarize(enriched)
    img = draw_top_compounds(enriched.sort_values("standard_value"))
    print(f"Saved top-compound grid image to {img}")


if __name__ == "__main__":
    main()
