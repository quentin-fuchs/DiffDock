"""
Audit PoseBusters benchmark files to identify format issues that could explain
the 29% top-1 vs 50% paper target.

Checks per complex:
  1. Atom count match: _ligand.sdf == _ligands.sdf (mol 0) == _ligand_start_conf.sdf
  2. Ligand centroid vs protein Cα centroid distance (flag if >30Å = displaced)
  3. SMILES match between _ligand.sdf and _ligand_start_conf.sdf
  4. Whether _ligands.sdf mol 0 atoms == _ligand.sdf atoms (alignment sanity)
  5. Protein Cα count and chain count
  6. Counts of inf RMSD in existing results and which complexes they are

Also compares our RMSD results (merged npy) vs what we expect and flags
complexes where top-1 RMSD is inf (spyrmsd timeout or error).
"""

import os, sys, warnings
import numpy as np
from collections import defaultdict

warnings.filterwarnings('ignore')
from rdkit import Chem
from rdkit.Chem import AllChem, RemoveAllHs

DATA_DIR   = "/home/qf226/MProject/DiffDock/data/posebusters_benchmark_set"
MERGED_DIR = "/home/qf226/rds/hpc-work/results/pb_evaluate_merged"
OUT_PATH   = "/home/qf226/MProject/DiffDock/analysis/pb_file_audit_results.txt"

# ── helpers ──────────────────────────────────────────────────────────────────

def mol_from_sdf(path, mol_idx=0):
    supp = Chem.SDMolSupplier(path, sanitize=False, removeHs=False)
    mols = list(supp)
    if mol_idx >= len(mols) or mols[mol_idx] is None:
        return None
    mol = mols[mol_idx]
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    return RemoveAllHs(mol)

def mol_count(path):
    try:
        return sum(1 for _ in Chem.SDMolSupplier(path, sanitize=False, removeHs=False))
    except Exception:
        return -1

def mol_smiles(mol):
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None

def ca_positions(pdb_path):
    """Return Cα atom coordinates as array (N, 3)."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[12:16].strip() == "CA":
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    pass
    return np.array(coords) if coords else None

def ligand_centroid(mol):
    if mol is None or mol.GetNumConformers() == 0:
        return None
    return mol.GetConformer().GetPositions().mean(axis=0)

# ── load existing RMSD results ────────────────────────────────────────────────

print("Loading existing RMSD results...")
rmsds   = np.load(os.path.join(MERGED_DIR, "rmsds.npy"))       # (305, 40)
names   = np.load(os.path.join(MERGED_DIR, "complex_names.npy"), allow_pickle=True)
confs   = np.load(os.path.join(MERGED_DIR, "confidences.npy")) # (305, 40)

conf_top1 = np.argmax(confs, axis=1)
top1_rmsd = rmsds[np.arange(len(rmsds)), conf_top1]
oracle    = np.min(rmsds, axis=1)

# complexes with inf in ANY sample
has_inf = np.any(~np.isfinite(rmsds), axis=1)
inf_oracle = ~np.isfinite(oracle)
inf_top1   = ~np.isfinite(top1_rmsd)
eval_names_set = set(names)

print(f"  Total complexes evaluated: {len(names)}")
print(f"  Complexes with ≥1 inf RMSD sample: {has_inf.sum()}")
print(f"  Oracle is inf: {inf_oracle.sum()}")
print(f"  Top-1 is inf: {inf_top1.sum()}")
print(f"  Top-1 < 2Å: {(top1_rmsd < 2).sum()} / {len(top1_rmsd)}  ({100*(top1_rmsd<2).mean():.2f}%)")
print(f"  Oracle < 2Å: {(oracle < 2).sum()} / {len(oracle)}  ({100*(oracle<2).mean():.2f}%)")

# ── per-complex file audit ────────────────────────────────────────────────────

complexes = sorted(os.listdir(DATA_DIR))
print(f"\nAuditing {len(complexes)} complexes in {DATA_DIR}...\n")

issues = defaultdict(list)  # name → list of issue strings
stats = defaultdict(int)

atom_count_mismatches = []
displaced = []
smiles_mismatches = []
inf_rmsd_complexes = []

for name in complexes:
    d = os.path.join(DATA_DIR, name)
    lig_path       = os.path.join(d, f"{name}_ligand.sdf")
    ligs_path      = os.path.join(d, f"{name}_ligands.sdf")
    start_path     = os.path.join(d, f"{name}_ligand_start_conf.sdf")
    prot_path      = os.path.join(d, f"{name}_protein.pdb")

    # --- file existence ---
    for p, label in [(lig_path, "ligand.sdf"), (ligs_path, "ligands.sdf"),
                     (start_path, "ligand_start_conf.sdf"), (prot_path, "protein.pdb")]:
        if not os.path.exists(p):
            issues[name].append(f"MISSING {label}")
            stats["missing_file"] += 1

    if not all(os.path.exists(p) for p in [lig_path, ligs_path, start_path, prot_path]):
        continue

    # --- molecule loading ---
    mol_lig    = mol_from_sdf(lig_path)
    mol_ligs0  = mol_from_sdf(ligs_path, 0)
    mol_start  = mol_from_sdf(start_path)
    n_copies   = mol_count(ligs_path)

    failed = []
    if mol_lig is None:   failed.append("ligand.sdf parse")
    if mol_ligs0 is None: failed.append("ligands.sdf[0] parse")
    if mol_start is None: failed.append("ligand_start_conf.sdf parse")
    if failed:
        for f in failed:
            issues[name].append(f"PARSE FAIL {f}")
            stats["parse_fail"] += 1
        continue

    # --- atom count checks ---
    na_lig   = mol_lig.GetNumAtoms()
    na_ligs0 = mol_ligs0.GetNumAtoms()
    na_start = mol_start.GetNumAtoms()

    if na_lig != na_ligs0:
        issues[name].append(f"ATOM COUNT ligand({na_lig}) != ligands[0]({na_ligs0})")
        stats["atom_count_mismatch"] += 1
        atom_count_mismatches.append((name, na_lig, na_ligs0, na_start))

    if na_lig != na_start:
        issues[name].append(f"ATOM COUNT ligand({na_lig}) != start_conf({na_start})")
        stats["start_conf_mismatch"] += 1
        if (name, na_lig, na_ligs0, na_start) not in atom_count_mismatches:
            atom_count_mismatches.append((name, na_lig, na_ligs0, na_start))

    # --- SMILES check: ligand.sdf vs ligand_start_conf.sdf ---
    smi_lig   = mol_smiles(mol_lig)
    smi_start = mol_smiles(mol_start)
    if smi_lig and smi_start and smi_lig != smi_start:
        issues[name].append(f"SMILES MISMATCH ligand vs start_conf")
        stats["smiles_mismatch"] += 1
        smiles_mismatches.append((name, smi_lig, smi_start))

    # --- displaced ligand check ---
    ca = ca_positions(prot_path)
    if ca is None or len(ca) == 0:
        issues[name].append("NO CA ATOMS in protein.pdb")
        stats["no_ca"] += 1
    else:
        prot_center = ca.mean(axis=0)
        lig_cen = ligand_centroid(mol_lig)
        if lig_cen is not None:
            dist = np.linalg.norm(lig_cen - prot_center)
            if dist > 30.0:
                issues[name].append(f"DISPLACED dist={dist:.1f}Å")
                stats["displaced"] += 1
                displaced.append((name, dist))

    # --- check if this complex has inf RMSD in results ---
    if name in eval_names_set:
        idx = np.where(names == name)[0]
        if len(idx) > 0:
            i = idx[0]
            n_inf = (~np.isfinite(rmsds[i])).sum()
            if n_inf > 0:
                inf_rmsd_complexes.append((name, n_inf, top1_rmsd[i], oracle[i], n_copies))

    stats["ok"] += 1

# ── print summary ─────────────────────────────────────────────────────────────

lines = []
lines.append("=" * 70)
lines.append("PoseBusters file audit summary")
lines.append("=" * 70)
lines.append(f"Total complexes in benchmark dir: {len(complexes)}")
lines.append(f"  Missing files:        {stats['missing_file']}")
lines.append(f"  Parse failures:       {stats['parse_fail']}")
lines.append(f"  Atom count mismatch (ligand.sdf vs ligands.sdf[0]):  {stats['atom_count_mismatch']}")
lines.append(f"  Atom count mismatch (ligand.sdf vs start_conf.sdf):  {stats['start_conf_mismatch']}")
lines.append(f"  SMILES mismatch (ligand vs start_conf):  {stats['smiles_mismatch']}")
lines.append(f"  Displaced ligands (centroid >30Å from Cα center):    {stats['displaced']}")
lines.append(f"  No Cα atoms in protein:  {stats['no_ca']}")
lines.append("")

if atom_count_mismatches:
    lines.append("--- Atom count mismatches ---")
    for name, na_lig, na_ligs0, na_start in atom_count_mismatches[:30]:
        lines.append(f"  {name}  ligand={na_lig}  ligands[0]={na_ligs0}  start_conf={na_start}")
    lines.append("")

if smiles_mismatches:
    lines.append("--- SMILES mismatches (ligand.sdf vs start_conf.sdf) ---")
    for name, s1, s2 in smiles_mismatches[:20]:
        lines.append(f"  {name}")
        lines.append(f"    ligand:     {s1[:80]}")
        lines.append(f"    start_conf: {s2[:80]}")
    lines.append("")

if displaced:
    lines.append("--- Displaced ligands (>30Å from protein Cα centroid) ---")
    for name, dist in sorted(displaced, key=lambda x: -x[1])[:30]:
        in_eval = "in_eval" if name in eval_names_set else "NOT_in_eval"
        lines.append(f"  {name}  dist={dist:.1f}Å  [{in_eval}]")
    lines.append("")

if inf_rmsd_complexes:
    lines.append("--- Complexes with inf RMSD in evaluation results ---")
    for name, n_inf, t1, orc, n_copies in sorted(inf_rmsd_complexes, key=lambda x: -x[1])[:40]:
        lines.append(f"  {name}  n_inf_samples={n_inf}/40  top1={t1:.2f}  oracle={orc:.2f}  n_copies={n_copies}")
    lines.append("")

lines.append("--- All per-complex issues ---")
for name, islist in sorted(issues.items()):
    lines.append(f"  {name}: {'; '.join(islist)}")

output = "\n".join(lines)
print(output)
with open(OUT_PATH, "w") as f:
    f.write(output)
print(f"\nResults saved to {OUT_PATH}")
