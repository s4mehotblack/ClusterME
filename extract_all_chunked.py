#!/usr/bin/env python3
import pandas as pd
import os
import sys
import argparse

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Extract GTEx expression data for genes of interest")
    
    parser.add_argument("--genes-file", required=True, help="File containing list of genes to extract")
    
    parser.add_argument("--sample-attributes", 
                        default="GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                        help="GTEx sample attributes file (default: GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt)")
    
    parser.add_argument("--gtex-tpm-file", 
                        default="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz",
                        help="GTEx TPM matrix file (default: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz)")
    
    parser.add_argument("--output-file", 
                        default="gtex_expression_all_tissues.csv",
                        help="Output CSV file (default: gtex_expression_all_tissues.csv)")
    
    parser.add_argument("--dry-run", 
                        action="store_true",
                        help="Dry run: report statistics without processing")

    args = parser.parse_args()

    try:
        # ──────────────────────────────────────────
        # LOAD GENE LIST
        # ──────────────────────────────────────────
        if not os.path.exists(args.genes_file):
            print(f"❌ Error: Genes file not found: {args.genes_file}")
            sys.exit(1)

        with open(args.genes_file) as f:
            genes_of_interest = {line.strip() for line in f if line.strip()}

        print(f"Loaded {len(genes_of_interest)} genes")

        # ──────────────────────────────────────────
        # LOAD SAMPLE ATTRIBUTES (GTEx)
        # ──────────────────────────────────────────
        if not os.path.exists(args.sample_attributes):
            print(f"❌ Error: Sample attributes file not found: {args.sample_attributes}")
            sys.exit(1)

        attrs = pd.read_csv(
            args.sample_attributes,
            sep="\t",
            usecols=["SAMPID", "SMTSD"]
        )

        # Map: tissue → sample IDs
        tissue_to_samples = (
            attrs.groupby("SMTSD")["SAMPID"]
            .apply(list)
            .to_dict()
        )

        print(f"Found {len(tissue_to_samples)} tissue types")

        # ──────────────────────────────────────────
        # DRY RUN MODE
        # ──────────────────────────────────────────
        if args.dry_run:
            print(f"\nDRY RUN SUMMARY:")
            print(f"  Genes to process: {len(genes_of_interest)}")
            print(f"  Tissue types found: {len(tissue_to_samples)}")
            print(f"  GTEx TPM file: {args.gtex_tpm_file}")
            print(f"  Output would be written to: {args.output_file}")
            print(f"  Chunksize: 2000 genes per chunk")
            print("\nDry run complete. No files were modified.")
            return

        # ──────────────────────────────────────────
        # OUTPUT FILE
        # ──────────────────────────────────────────
        out_file = args.output_file

        if os.path.exists(out_file):
            try:
                os.remove(out_file)
            except OSError as e:
                print(f"❌ Error: Could not remove existing output file {out_file}: {e}")
                sys.exit(1)

        header_written = False

        # ──────────────────────────────────────────
        # STREAM GTEx TPM MATRIX
        # ──────────────────────────────────────────
        if not os.path.exists(args.gtex_tpm_file):
            print(f"❌ Error: GTEx TPM file not found: {args.gtex_tpm_file}")
            sys.exit(1)

        chunksize = 2000

        reader = pd.read_csv(
            args.gtex_tpm_file,
            sep="\t",
            skiprows=2,
            chunksize=chunksize,
            low_memory=True
        )

        for i, chunk in enumerate(reader):
            print(f"Processing chunk {i}")

            # Standardize gene column
            chunk.rename(columns={"Description": "gene"}, inplace=True)

            # Filter genes
            sub = chunk[chunk["gene"].isin(genes_of_interest)]
            if sub.empty:
                continue

            sub = sub.copy()

            rows = []

            for tissue, samples in tissue_to_samples.items():
                common = [s for s in samples if s in sub.columns]
                if not common:
                    continue

                mean_tpm = sub[common].mean(axis=1)
                median_tpm = sub[common].median(axis=1)
                n_samples = len(common)

                tmp = pd.DataFrame({
                    "gene": sub["gene"],
                    "tissue": tissue,
                    "mean_tpm": mean_tpm,
                    "median_tpm": median_tpm,
                    "n_samples": n_samples
                })

                rows.append(tmp)

            if not rows:
                continue

            out = pd.concat(rows, ignore_index=True)

            try:
                out.to_csv(
                    out_file,
                    index=False,
                    mode="a",
                    header=not header_written
                )
                header_written = True
            except Exception as e:
                print(f"❌ Error writing to output file {out_file}: {e}")
                sys.exit(1)

        print(f"\nDone. Output written to: {out_file}\n")

    except KeyboardInterrupt:
        print("\n🛑 Extraction interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
