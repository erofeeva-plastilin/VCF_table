import pandas as pd
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

primers = pd.read_csv("primers.tsv", sep="\t")
primers.set_index("region", inplace=True)


amplicon_seqs = {
    region: row.tolist()
    for region, row in primers.iterrows()
}


with open("list_fq.tsv") as f:
    fq_files = sorted([line.strip() for line in f if line.strip()])

from collections import defaultdict

fq_pairs = defaultdict(dict)
for fq in fq_files:
    name = Path(fq).name
    sample = name.rsplit("_R", 1)[0]
    if "_R1" in fq:
        fq_pairs[sample]["R1"] = fq
    elif "_R2" in fq:
        fq_pairs[sample]["R2"] = fq

samples = sorted(fq_pairs.keys())
amplicons = sorted(amplicon_seqs.keys())

def count_reads(sample, r1, r2, region, seqs):
    pattern = "|".join(seqs)
    try:
        r1_count = int(subprocess.check_output(
            f"zcat {r1} | grep -E -c '{pattern}'", shell=True, text=True
        ).strip())
        r2_count = int(subprocess.check_output(
            f"zcat {r2} | grep -E -c '{pattern}'", shell=True, text=True
        ).strip())
        return (region, sample, r1_count + r2_count)
    except subprocess.CalledProcessError:
        return (region, sample, 0)

MAX_THREADS = 12
results = []

with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
    futures = []
    for sample, pair in fq_pairs.items():
        r1, r2 = pair.get("R1"), pair.get("R2")
        if not r1 or not r2:
            print(f"[WARN] Пропущен неполный сэмпл: {sample}")
            continue
        for region in amplicons:
            seqs = amplicon_seqs[region]
            futures.append(executor.submit(count_reads, sample, r1, r2, region, seqs))

    for future in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
        results.append(future.result())

df = pd.DataFrame(results, columns=["amplicon", "sample", "count"])
pivot = df.pivot(index="amplicon", columns="sample", values="count").fillna(0).astype(int)
pivot.to_csv("amplicon_fastq_counts.tsv", sep="\t")

print("amplicon_fastq_counts.tsv")
