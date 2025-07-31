import sys
import os

def process_vcf(input_file):
    try: 
        # The .vcf extension is removed from the input file name to obtain a base name
        # Using this base name, two output file names are formed: one for genotypes and one for allele depths
        base_name = os.path.splitext(input_file)[0]
        gt_output = f"{base_name}_GT.tsv"
        ad_output = f"{base_name}_AD.tsv"
        
        # The input VCF file is opened and read line by line
        # If a line starts with #, it is considered a header line
        # If the line contains sample names (starts with #CHROM), 
        # the column and sample names are extracted to form the headers for the output tables
        # A flag is declared at the beginning to track whether the allele depth (AD) file has been created
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
                    
                    # Splitting the contents of the line
                    fields = line.strip().split('\t')
                    chrom, pos, _, ref, alt, qual, _, _, fmt = fields[:9]
                    sample_data = fields[9:]

                    # The indices of the GT and AD fields in the format are defined
                    fmt_fields = fmt.split(':')
                    
                    # If the GT field is present, its index is stored in the variable gt_idx, otherwise the value will be -1
                    gt_idx = fmt_fields.index("GT") if "GT" in fmt_fields else -1
                    
                    # Similarly for the AD field. This allows the correct positions to be used later 
                    # to extract GT and AD values from the sample strings
                    ad_idx = fmt_fields.index("AD") if "AD" in fmt_fields else -1

                    # If the AD field is detected for the first time, an AD file is created and a header is written
                    if ad_idx != -1 and not has_ad:
                        has_ad = True
                        ad_header = ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL'] + samples
                        ad_file = open(ad_output, 'w')
                        ad_file.write('\t'.join(ad_header) + '\n')

                    # For each row of data from the VCF file, two arrays are created:
                    # one for genotype data (GT) and one for allele depth (AD), if present (if the has_ad flag is set to True)
                    gt_row = [chrom, pos, ref, alt, qual]
                    ad_row = [chrom, pos, ref, alt, qual] if has_ad else None

                    # Data extraction for each sample
                    # Sample data fields are separated by the : symbol
                    # For example, if the sample data is in the format 0/1:10,5,
                    # after splitting you will get a list: ['0/1', '10,5']
                    for sample in sample_data:
                        sample_fields = sample.split(':')
                        
                        # Extracting and adding genotype (GT) to a string
                        gt = sample_fields[gt_idx] if gt_idx != -1 else "."
                        gt_row.append(gt)

                        # Extracting and adding allelic depth (AD) to a string
                        if has_ad:
                            if ad_idx != -1 and ad_idx < len(sample_fields):
                                ad = sample_fields[ad_idx]
                                if ad != ".":
                                    ad_values = ad.split(',')
                                    ad = "/".join(ad_values)
                                else:
                                    ad = "./."
                            else:
                                ad = "./."
                            ad_row.append(ad)

                    # Writing data to output files
                    gt_file.write('\t'.join(gt_row) + '\n')
                    if has_ad:
                        ad_file.write('\t'.join(ad_row) + '\n')

        # Closing the AD file if it was created
        if has_ad and ad_file:
            ad_file.close()
    except Exception as e:
        print(f"Error during file processing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Checking the number of command line arguments so that there is no error when multiple files are submitted
    if len(sys.argv) != 2:
        sys.exit(1)
    process_vcf(sys.argv[1])
