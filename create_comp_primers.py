import pandas as pd

def complement(seq):
    comp_dict = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(comp_dict)

def reverse(seq):
    return seq[::-1]

def reverse_complement(seq):
    return complement(reverse(seq))

input_file = "all.txt"
output_file = "primers_with_variants1.tsv"

df = pd.read_csv(input_file, sep=r"\s+", engine="python")

for col in ["left", "right"]:
    df[f"{col}_complement"] = df[col].apply(complement)
    df[f"{col}_reverse"] = df[col].apply(reverse)
    df[f"{col}_revcomp"] = df[col].apply(reverse_complement)

df.to_csv(output_file, sep="\t", index=False)

print(f"Готово! Файл сохранён как {output_file}")
