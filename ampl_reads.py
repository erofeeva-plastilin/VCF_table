import subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

with open("cram_list.tsv") as f:
    cram_files = [line.strip() for line in f if line.strip()]

with open("ampl_list.tsv") as f:
    amplicons = [line.strip() for line in f if line.strip()]

def get_cram_references(cram_path):
    refs = []
    try:
        header = subprocess.check_output(["samtools", "view", "-H", cram_path], text=True)
        for line in header.splitlines():
            if line.startswith("@SQ"):
                for field in line.split("\t"):
                    if field.startswith("SN:"):
                        refs.append(field.replace("SN:", ""))
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to read CRAM header from {cram_path}: {e}")
    return refs

cram_refs = {cram: set(get_cram_references(cram)) for cram in cram_files}

def count_reads(cram, amp):
    if amp not in cram_refs[cram]:
        return (amp, cram, 0)
    try:
        count = subprocess.check_output(
            ["samtools", "view", "-c", cram, amp],
            stderr=subprocess.DEVNULL,
            text=True
        ).strip()
        return (amp, cram, int(count))
    except subprocess.CalledProcessError:
        return (amp, cram, 0)

results_dict = {amp: {} for amp in amplicons}
MAX_THREADS = 20

with ThreadPoolExecutor(max_workers=MAX_THREADS) as executor:
    futures = {
        executor.submit(count_reads, cram, amp): (cram, amp)
        for cram in cram_files
        for amp in amplicons
    }

    for future in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
        amp, cram, count = future.result()
        results_dict[amp][cram] = count

df = pd.DataFrame.from_dict(results_dict, orient="index", columns=cram_files)
df.index.name = "amplicon"
df = df[[c for c in cram_files]]
df.columns = [c.split("/")[-1] for c in df.columns]

df.to_csv("amplicon_read_counts.tsv", sep="\t")
print("amplicon_read_counts.tsv")
