# Parse the Carrot-Zhang 2020 ancestry table into a tidy per-patient parquet.
# The table is a Cell supplement (Table S1). Download it, find the header row
# under the title banner, and pull the consensus ancestry call plus the
# admixture percentages.

import os
import requests
import pandas as pd

DATA_ROOT = os.environ["DATA_ROOT"]
raw = os.path.join(DATA_ROOT, "raw")
interim = os.path.join(DATA_ROOT, "interim")
os.makedirs(raw, exist_ok=True)

# the supplement lives on the Cell/ScienceDirect CDN under the article PII
pii = "S1535610820302117"
url = "https://ars.els-cdn.com/content/image/1-s2.0-" + pii + "-mmc2.xlsx"
fn = os.path.join(raw, "carrotzhang_mmc2.xlsx")
if not os.path.exists(fn):
    r = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=90)
    r.raise_for_status()
    open(fn, "wb").write(r.content)

# the per-patient calls are on the "S1 Calls per Patient" sheet, with the real
# header a row or two below the title banner
xl = pd.ExcelFile(fn)
sheet = "S1 Calls per Patient"
probe = xl.parse(sheet, header=None, nrows=6)
header_row = None
for i in range(len(probe)):
    cells = [str(v).lower() for v in probe.iloc[i].tolist()]
    if any("consensus_ancestry" in c for c in cells):
        header_row = i
        break
if header_row is None:
    header_row = 1

df = xl.parse(sheet, header=header_row)
df.columns = [str(c).strip() for c in df.columns]
df["patient12"] = df["patient"].astype(str).str.slice(0, 12)

# normalize the admixture percentage column names
rename = {}
for c in df.columns:
    cl = c.lower()
    if "admixture" in cl and "afr" in cl:
        rename[c] = "pct_afr"
    elif "admixture" in cl and "eur" in cl:
        rename[c] = "pct_eur"
df = df.rename(columns=rename)
df = df.rename(columns={"patient": "patient_barcode"})

keep = ["patient_barcode", "patient12", "consensus_ancestry", "tumor_type"]
for c in ["pct_afr", "pct_eur"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
        keep.append(c)

out = df[keep].copy()
out.to_parquet(os.path.join(interim, "tcga_ancestry.parquet"), index=False)
out.to_csv(os.path.join(interim, "tcga_ancestry.csv"), index=False)
print("saved ancestry table:", len(out), "patients")
