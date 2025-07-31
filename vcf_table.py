import sys
import os

def process_vcf(input_file):
    try: 
        base_name = os.path.splitext(input_file)[0]
        gt_output = f"{base_name}_GT.tsv"
        ad_output = f"{base_name}_AD.tsv"
        
        with open(input_file, 'r') as vcf:
            has_ad = False
            with open(gt_output, 'w') as gt_file:
                ad_file = None
                for line in vcf:
                    if line.startswith('#'):
                        if line.startswith("#CHROM"):
                            header = line.strip().split('\t')
                            samples = header[9:]
                            gt_header = ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL'] + samples
                            gt_file.write('\t'.join(gt_header) + '\n')
                        continue
                    
                    fields = line.strip().split('\t')
                    chrom, pos, _, ref, alt, qual, _, _, fmt = fields[:9]
                    sample_data = fields[9:]

                    fmt_fields = fmt.split(':')
                    fmt_dict = {key: i for i, key in enumerate(fmt_fields)}
                    gt_idx = fmt_dict.get("GT", -1)
                    ad_idx = fmt_dict.get("AD", -1)

                    if ad_idx != -1 and not has_ad:
                        has_ad = True
                        ad_header = ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL'] + samples
                        ad_file = open(ad_output, 'w')
                        ad_file.write('\t'.join(ad_header) + '\n')

                    gt_row = [chrom, pos, ref, alt, qual]
                    ad_row = [chrom, pos, ref, alt, qual] if has_ad else None

                    for sample in sample_data:
                        sample_fields = sample.split(':')

                        # GT
                        gt = sample_fields[gt_idx] if gt_idx != -1 and gt_idx < len(sample_fields) else "."
                        gt_row.append(gt)

                        # AD
                        if has_ad:
                            ad = "./."
                            if ad_idx != -1 and ad_idx < len(sample_fields):
                                ad_raw = sample_fields[ad_idx]
                                if ad_raw != ".":
                                    ad_values = ad_raw.split(',')
                                    ad = "/".join(ad_values)
                            ad_row.append(ad)

                    gt_file.write('\t'.join(gt_row) + '\n')
                    if has_ad:
                        ad_file.write('\t'.join(ad_row) + '\n')

        if has_ad and ad_file:
            ad_file.close()
    except Exception as e:
        print(f"Error during file processing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.vcf")
        sys.exit(1)
    process_vcf(sys.argv[1])
