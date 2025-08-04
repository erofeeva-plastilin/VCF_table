import pandas as pd
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

primers = pd.read_csv("primers.tsv", sep="\t")
primers.set_index("region", inplace=True)
amplicon_seqs = {region: row.tolist() for region, row in primers.iterrows()}

with open("cram_list.tsv") as f:
    cram_files = sorted([line.strip() for line in f if line.strip()])

sample_crams = {Path(f).stem.replace(".cram", ""): f for f in cram_files}
samples = sorted(sample_crams.keys())
amplicons = sorted(amplicon_seqs.keys())

def count_reads_start_match(sample, cram_path, region, seqs):
    pattern = "^(" + "|".join(seqs) + ")"
    try:
        cmd = f"samtools view -f 4 {cram_path} | grep -E -c '{pattern}'"
        count = int(subprocess.check_output(cmd, shell=True, text=True).strip())
        return (region, sample, count)
    except subprocess.CalledProcessError:
        return (region, sample, 0)

MAX_THREADS = 20
results = []

with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
    futures = []
    for sample, cram in sample_crams.items():
        for region in amplicons:
            seqs = amplicon_seqs[region]
            futures.append(executor.submit(count_reads_start_match, sample, cram, region, seqs))

    for future in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
        results.append(future.result())

df = pd.DataFrame(results, columns=["amplicon", "sample", "count"])
pivot = df.pivot(index="amplicon", columns="sample", values="count").fillna(0).astype(int)
pivot.to_csv("start_amplicon_cram_startmatch_counts.tsv", sep="\t")

print("start_amplicon_cram_startmatch_counts.tsv")
