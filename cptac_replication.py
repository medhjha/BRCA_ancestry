# CPTAC replication and protein corroboration.
# For each clean-core edge, recompute the EUR-minus-AFR correlation difference in
# CPTAC and compare direction and magnitude to TCGA. We use the raw
# EUR-minus-AFR correlation difference in both cohorts so the comparison does not
# depend on the sign convention of any single statistic.

import os
import cptac
import pandas as pd
import numpy as np
from scipy import stats

RESULTS = os.environ.get("RESULTS", "results")

# clean-core edges exported from R (columns G1, G2, EUR_cor, AFR_cor, zScoreDiff)
core = pd.read_csv(os.path.join(RESULTS, "clean_core_edges_for_replication.csv"))
core["tcga_delta"] = core["EUR_cor"] - core["AFR_cor"]

# CPTAC RNA and its own consensus ancestry calls
br = cptac.Brca()
rna = br.get_dataframe("transcriptomics", "bcm")
if isinstance(rna.columns, pd.MultiIndex):
    rna.columns = rna.columns.get_level_values(0)
rna = rna.loc[:, ~rna.columns.duplicated()]

anc_file = os.path.join(os.path.dirname(cptac.__file__),
                        "data/harmonized-all_cancers/washu_mssm_consensus_ancestries.tsv.gz")
anc = pd.read_csv(anc_file, sep="\t", compression="gzip")
anc = anc[anc["cancer_type"] == "BRCA"].set_index("case_id")["consensus_pred_ancestry"]
eur = [s for s in rna.index if anc.get(s) == "EUR"]
afr = [s for s in rna.index if anc.get(s) == "AFR"]
print("CPTAC RNA: EUR", len(eur), "AFR", len(afr))

def cptac_delta(g1, g2, mat, e, a):
    if g1 not in mat.columns or g2 not in mat.columns:
        return np.nan
    ee = mat.loc[e, [g1, g2]].dropna()
    aa = mat.loc[a, [g1, g2]].dropna()
    if len(ee) < 10 or len(aa) < 8:
        return np.nan
    return ee[g1].corr(ee[g2]) - aa[g1].corr(aa[g2])

core["cpd"] = [cptac_delta(r.G1, r.G2, rna, eur, afr) for r in core.itertuples()]
tested = core.dropna(subset=["cpd"]).copy()

def report(df, label):
    n = len(df)
    same = (np.sign(df["tcga_delta"]) == np.sign(df["cpd"]))
    conc = same.mean()
    p = stats.binomtest(int(same.sum()), n, 0.5, alternative="greater").pvalue
    r, rp = stats.spearmanr(df["tcga_delta"], df["cpd"])
    print(f"{label}: n={n}, {100*conc:.0f}% concordant (p={p:.1e}), Spearman={r:.2f}")

report(tested.head(200), "RNA top-200")
report(tested, "RNA all")
tested.to_csv(os.path.join(RESULTS, "cptac_replication.csv"), index=False)

# protein corroboration on the same patients
prot = br.get_dataframe("proteomics", "bcm")
if isinstance(prot.columns, pd.MultiIndex):
    prot.columns = prot.columns.get_level_values(0)
prot = prot.loc[:, ~prot.columns.duplicated()]
pe = [s for s in prot.index if anc.get(s) == "EUR"]
pa = [s for s in prot.index if anc.get(s) == "AFR"]
core["prot_delta"] = [cptac_delta(r.G1, r.G2, prot, pe, pa) for r in core.itertuples()]
pt = core.dropna(subset=["prot_delta"]).copy()
same = (np.sign(pt["tcga_delta"]) == np.sign(pt["prot_delta"]))
r, _ = stats.spearmanr(pt["tcga_delta"], pt["prot_delta"])
print(f"protein: n={len(pt)}, {100*same.mean():.0f}% concordant, Spearman={r:.2f}")
pt.to_csv(os.path.join(RESULTS, "cptac_protein_corroboration.csv"), index=False)
