"""
Find genes from a master list that exist in multiple gene set/cluster files.
Usage: python gene_intersection.py -m master_list.txt -s gene_set1.txt gene_set2.txt ...
"""

import sys
import argparse
from difflib import SequenceMatcher

def read_genes(filename):
    """Read genes from a file, one per line, removing whitespace."""
    with open(filename, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def find_exact_matches(master_genes, gene_set):
    """Find exact matches between master list and gene set."""
    return master_genes.intersection(gene_set)

def find_prefix_matches(master_genes, gene_set, min_length=3):
    """Find genes where the beginning of the gene name matches."""
    matches = set()
    for master_gene in master_genes:
        for set_gene in gene_set:
            if len(master_gene) >= min_length and len(set_gene) >= min_length:
                if master_gene.startswith(set_gene) or set_gene.startswith(master_gene):
                    matches.add(f"{master_gene}~{set_gene}")
                    break
    return matches

def levenshtein_distance(s1, s2):
    """Calculate the Levenshtein distance between two strings."""
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    
    if len(s2) == 0:
        return len(s1)
    
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]

def find_edit_distance_matches(master_genes, gene_set, max_distance=2):
    """Find genes within a certain edit distance."""
    matches = set()
    for master_gene in master_genes:
        for set_gene in gene_set:
            distance = levenshtein_distance(master_gene.upper(), set_gene.upper())
            if distance <= max_distance:
                matches.add(f"{master_gene}~{set_gene}(d={distance})")
                break
    return matches

def find_substring_matches(master_genes, gene_set, min_length=4):
    """Find genes that contain each other as substrings."""
    matches = set()
    for master_gene in master_genes:
        for set_gene in gene_set:
            master_upper = master_gene.upper()
            set_upper = set_gene.upper()
            if len(master_gene) >= min_length and len(set_gene) >= min_length:
                if master_upper in set_upper or set_upper in master_upper:
                    matches.add(f"{master_gene}~{set_gene}")
                    break
    return matches

def find_similarity_matches(master_genes, gene_set, threshold=0.8):
    """Find genes with similarity ratio above threshold using SequenceMatcher."""
    matches = set()
    for master_gene in master_genes:
        for set_gene in gene_set:
            ratio = SequenceMatcher(None, master_gene.upper(), set_gene.upper()).ratio()
            if ratio >= threshold:
                matches.add(f"{master_gene}~{set_gene}(s={ratio:.2f})")
                break
    return matches

def find_ngram_matches(master_genes, gene_set, n=3, threshold=0.5):
    """Find genes with similar n-grams."""
    def get_ngrams(string, n):
        string = string.upper()
        return set(string[i:i+n] for i in range(len(string)-n+1))
    
    matches = set()
    for master_gene in master_genes:
        master_ngrams = get_ngrams(master_gene, n)
        if not master_ngrams:
            continue
            
        for set_gene in gene_set:
            set_ngrams = get_ngrams(set_gene, n)
            if not set_ngrams:
                continue
                
            # Jaccard similarity
            intersection = len(master_ngrams & set_ngrams)
            union = len(master_ngrams | set_ngrams)
            similarity = intersection / union if union > 0 else 0
            
            if similarity >= threshold:
                matches.add(f"{master_gene}~{set_gene}(ng={similarity:.2f})")
                break
    return matches

def main():
    parser = argparse.ArgumentParser(
        description='Find genes from a master list that exist in gene set/cluster files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Matching modes:
  exact      : Exact string match (default)
  prefix     : Matches if one gene name starts with the other
  substring  : Matches if one gene name contains the other
  edit       : Matches within specified edit distance (typos, insertions, deletions)
  similarity : Matches with similarity ratio above threshold
  ngram      : Matches based on shared character n-grams
        """
    )
    parser.add_argument('-m', '--master', required=True,
                        help='Master gene list file (one gene per line)')
    parser.add_argument('-s', '--sets', nargs='+', required=True,
                        help='Gene set/cluster files to check against master list')
    parser.add_argument('--mode', choices=['exact', 'prefix', 'substring', 'edit', 'similarity', 'ngram'],
                        default='exact',
                        help='Matching mode (default: exact)')
    
    # Mode-specific parameters
    parser.add_argument('--min-length', type=int, default=3,
                        help='Minimum length for prefix/substring matching (default: 3)')
    parser.add_argument('--max-distance', type=int, default=2,
                        help='Maximum edit distance for edit mode (default: 2)')
    parser.add_argument('--threshold', type=float, default=0.8,
                        help='Similarity threshold for similarity/ngram modes (default: 0.8)')
    parser.add_argument('--ngram-size', type=int, default=3,
                        help='N-gram size for ngram mode (default: 3)')
    
    args = parser.parse_args()
    
    master_file = args.master
    gene_set_files = args.sets
    
    # Read master list
    if not os.path.exists(master_file):
        print(f"❌ Error: Master file not found: {master_file}")
        sys.exit(1)

    try:
        master_genes = read_genes(master_file)
    except Exception as e:
        print(f"❌ Error reading master file: {e}")
        sys.exit(1)

    print(f"Master list ({master_file}) contains {len(master_genes)} genes")
    print(f"Matching mode: {args.mode}\n")
    
    # Check each gene set file
    for gene_set_file in gene_set_files:
        if not os.path.exists(gene_set_file):
            print(f"⚠ Warning: Gene set file not found: {gene_set_file}")
            continue

        try:
            gene_set = read_genes(gene_set_file)
        except Exception as e:
            print(f"❌ Error reading gene set file {gene_set_file}: {e}")
            continue
        
        # Find matches based on mode
        if args.mode == 'exact':
            common_genes = find_exact_matches(master_genes, gene_set)
        elif args.mode == 'prefix':
            common_genes = find_prefix_matches(master_genes, gene_set, args.min_length)
        elif args.mode == 'substring':
            common_genes = find_substring_matches(master_genes, gene_set, args.min_length)
        elif args.mode == 'edit':
            common_genes = find_edit_distance_matches(master_genes, gene_set, args.max_distance)
        elif args.mode == 'similarity':
            common_genes = find_similarity_matches(master_genes, gene_set, args.threshold)
        elif args.mode == 'ngram':
            common_genes = find_ngram_matches(master_genes, gene_set, args.ngram_size, args.threshold)
        
        # Output results
        if common_genes:
            genes_list = ','.join(sorted(common_genes))
            print(f"{gene_set_file} ({len(common_genes)} of {len(gene_set)} genes match): {genes_list}")
        else:
            print(f"{gene_set_file} (0 of {len(gene_set)} genes match): No matches")

if __name__ == "__main__":
    main()
