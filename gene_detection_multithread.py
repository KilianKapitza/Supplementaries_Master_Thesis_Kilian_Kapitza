from Bio import SeqIO
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial
import numpy as np
import pandas as pd
import argparse 
import hashlib
import sys
import re

hash_length = 19
appending_length = 18

def hash_kmer(sequence):
    """Generate a unique integer hash value for a k-mer sequence using SHA256."""
    hash_object = hashlib.sha256(sequence.encode('utf-8'))
    return int(hash_object.hexdigest()[:8], 16)  # Truncate to 8 hex digits for consistency

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                 'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                 'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))


def generate_reference_hashes(reference_file, include_reverse_complement=True):
    """Generate hash values for reference regions with various mutation types."""
    try:
        df = pd.read_csv(reference_file, sep='\t')
    except Exception as e:
        raise ValueError(f"Error reading the reference file '{reference_file}': {e}")
    
    reference_hashes = {}  # Initialize dictionary to store hashes
    removed_count = 0
    
    for _, row in df.iterrows():
        upstream = row['Upstream'].strip()
        ncf = row['NCF1_ungleich_BC'].strip()
        downstream = row['Downstream'].strip()
        mutation_type = row['Type'].strip()
        
        # Trim sequences according to new format: 
        # Use last 14 bases of upstream and first 14 bases of downstream
        if len(upstream) > appending_length:
            trimmed_upstream = upstream[-appending_length:]
        else:
            trimmed_upstream = upstream
            
        if len(downstream) > appending_length:
            trimmed_downstream = downstream[:appending_length]
        else:
            trimmed_downstream = downstream
        
        if mutation_type == "SNP":
            # For SNP: Create exactly 15 k-mers of exactly 15bp each
            if len(trimmed_upstream) == appending_length and len(ncf) == 1 and len(trimmed_downstream) == appending_length:
                # Create k-mers with SNP at each possible position (0-14)
                for snp_pos in range(hash_length):
                    # Calculate how many bases to take from upstream and downstream
                    upstream_bases = snp_pos
                    downstream_bases = hash_length - snp_pos - 1  # -1 for the SNP itself
                    
                    # Extract the required portions of the sequences
                    if upstream_bases > 0:
                        upstream_part = trimmed_upstream[-upstream_bases:]
                    else:
                        upstream_part = ""
                        
                    if downstream_bases > 0:
                        downstream_part = trimmed_downstream[:downstream_bases]
                    else:
                        downstream_part = ""
                    
                    # Combine to form the 15-mer with SNP at position snp_pos
                    mer = upstream_part + ncf + downstream_part
                    
                    if len(mer) == hash_length:
                        reference_hashes[mer] = hash_kmer(mer)
                    else:
                        removed_count += 1
                        print(f"Skipped mer of length {len(mer)}: {mer}")
            else:
                removed_count += 1
        
        elif mutation_type == "Insertion":
            # For Insertion: Handle insertion with trimmed sequences
            ncf_length = len(ncf)
            if ncf_length > hash_length:
                # Make 15-mers only out of the ncf region
                for i in range(ncf_length - hash_length + 1):
                    mer = ncf[i:i + hash_length]
                    if len(mer) == hash_length:
                        reference_hashes[mer] = hash_kmer(mer)
                    else:
                        removed_count += 1
            elif ncf_length == hash_length:
                # Create a single 15-mer with insertion in the middle
                mer = ncf 
                if len(mer) == hash_length:
                    reference_hashes[mer] = hash_kmer(mer)
                else:
                    removed_count += 1
            elif ncf_length > 1 and ncf_length < hash_length:
                # Create multiple 15-mers with insertion at different positions
                remaining_length = hash_length - ncf_length
                
                # Calculate max positions to shift the insertion
                max_shift = min(appending_length, remaining_length) + 1
                
                for shift in range(max_shift):
                    # Calculate how many bases to take from upstream and downstream
                    upstream_bases = shift
                    downstream_bases = remaining_length - shift
                    
                    # Get parts of sequences
                    if upstream_bases > 0:
                        upstream_part = trimmed_upstream[-upstream_bases:]
                    else:
                        upstream_part = ""
                        
                    if downstream_bases > 0:
                        downstream_part = trimmed_downstream[:downstream_bases]
                    else:
                        downstream_part = ""
                    
                    # Create k-mer
                    mer = upstream_part + ncf + downstream_part
                    
                    if len(mer) == hash_length:
                        reference_hashes[mer] = hash_kmer(mer)
                    else:
                        removed_count += 1
            else:
                removed_count += 1
                
        elif mutation_type == "Deletion":
            # For Deletion: Combine trimmed upstream and downstream
            combined = trimmed_upstream + trimmed_downstream
            # Get all possible 15-mers from the combined sequence
            for i in range(max(0, len(combined) - hash_length + 1)):
                mer = combined[i:i+hash_length]
                if len(mer) == hash_length:
                    reference_hashes[mer] = hash_kmer(mer)
                else:
                    removed_count += 1
        else:
            removed_count += 1

    # After generating all reference hashes, add reverse complements if requested
    if include_reverse_complement:
        original_kmers = list(reference_hashes.keys())
        for kmer in original_kmers:
            rev_comp = reverse_complement(kmer)
            if rev_comp not in reference_hashes:  # Avoid duplicates
                reference_hashes[rev_comp] = hash_kmer(rev_comp)
        
        print(f"Added {len(reference_hashes) - len(original_kmers)} reverse complement k-mers")
    
    return reference_hashes, removed_count



# Separate Histogramm-Funktion
def plot_histogram(matches_per_read, threshold, output_file):
    """Plot histogram of matches per read and save to file."""
    plt.figure(figsize=(10, 6))
    cap_value = 2 * threshold
    modified_matches = [min(m, cap_value) for m in matches_per_read]
    exceeding_count = sum(1 for m in matches_per_read if m > cap_value)
    bins = range(0, cap_value)
    counts, _, _ = plt.hist(modified_matches, bins=bins, alpha=0.7, 
                           color='blue', edgecolor='black')
    
    if exceeding_count > 0:
        counts[-2] += exceeding_count
        
    plt.yscale('log')
    plt.title('Distribution of Reference Matches per Read (Log Scale)')
    plt.xlabel('Number of Reference Matches')
    plt.ylabel('Number of Reads (log scale)')
    plt.grid(True, alpha=0.3)
    plt.axvline(threshold, color='red', linestyle='dashed', linewidth=2, 
               label=f'Threshold ({threshold})')
        
    if exceeding_count > 0:
        plt.annotate(f'+{exceeding_count} reads > {cap_value}', 
                    xy=(cap_value, counts[-2]), 
                    xytext=(cap_value-5, counts[-2]*1.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=1.5))
        
    plt.legend()
    plt.tight_layout()
    
    # Statistiken ausgeben
    print(f"Match statistics:")
    print(f"  Average matches per read: {sum(matches_per_read)/len(matches_per_read):.2f}")
    print(f"  Maximum matches in a read: {max(matches_per_read)}")
    print(f"  Reads with >={threshold} matches: {sum(1 for m in matches_per_read if m >= threshold)}")
    print(f"  Reads with >{cap_value} matches: {exceeding_count}")
    print(f"  Total reads processed: {len(matches_per_read)}")
    
    plt.savefig(output_file)
    print(f"Histogram saved to {output_file}")

def process_read_batch(read_batch, reference_hash_set, hash_length, threshold):
    """Process a batch of reads in parallel"""
    results = []
    for record in read_batch:
        read_sequence = str(record.seq)
        # Generate all possible 11-mers and their hashes
        read_hashes = [hash_kmer(read_sequence[i:i + hash_length]) 
                      for i in range(len(read_sequence) - hash_length + 1)]
        
         
        unique_read_hashes = set(read_hashes)
        # Count matches with reference hashes (each reference k-mer counts only once)
        match_count = sum(1 for h in unique_read_hashes if h in reference_hash_set)
        
        
        # Store result
        is_real_gene = match_count >= threshold
        if is_real_gene:
            record.description += f" real_gene={int(is_real_gene)} matches={match_count}"
            results.append((record, match_count, True))
        else:
            results.append((None, match_count, False))
            
    return results

# Neue Funktion zum Speichern von Match-Statistiken hinzufügen
def save_match_statistics(matches_per_read, stats_file):
    """Save match statistics to a file."""
    with open(stats_file, 'w') as f:
        for count in matches_per_read:
            f.write(f"{count}\n")
    print(f"Match statistics saved to {stats_file}")


def classify_reads(fastq_file, reference_hashes, threshold=20, plot_histo=False, num_processes=None, stats_file=None):
    """Classify reads from the FASTQ file based on SNP regions using multiprocessing."""
    if num_processes is None:
        num_processes = mp.cpu_count()  # Use all available CPU cores by default
    
    print(f"Using {num_processes} processes for parallel processing")
    
    # Convert references to a set for faster lookup
    reference_hash_set = set(reference_hashes.values())
    
    # Handle input as stdin when fastq_file is '-'
    if fastq_file == '-':
        print("Reading from standard input...")
        # SeqIO kann direkt von sys.stdin lesen
        all_reads = list(SeqIO.parse(sys.stdin, "fastq"))
    else:
        all_reads = list(SeqIO.parse(fastq_file, "fastq"))
    
    total_reads = len(all_reads)
    print(f"Processing {total_reads} reads...")
    
    # Split reads into batches for parallel processing
    batch_size = max(100, total_reads // (num_processes * 10))  # Adjust batch size for efficiency
    read_batches = [all_reads[i:i + batch_size] for i in range(0, total_reads, batch_size)]
    
    # Create a partial function with fixed parameters
    process_func = partial(process_read_batch, 
                          reference_hash_set=reference_hash_set,
                          hash_length=hash_length,
                          threshold=threshold)
    
    # Process batches in parallel
    classified_reads = []
    matches_per_read = []
    
    with mp.Pool(processes=num_processes) as pool:
        # Process batches and collect results
        for batch_results in pool.imap_unordered(process_func, read_batches):
            for record, match_count, is_classified in batch_results:
                matches_per_read.append(match_count)
                if is_classified:
                    classified_reads.append(record)
    
    if stats_file:
        save_match_statistics(matches_per_read, stats_file)

    # Plot histogram if requested
    if plot_histo and args.plot:
        plot_histogram(matches_per_read, threshold, args.plot)

    return classified_reads

def extract_match_counts(fastq_file):
    """Extract match counts from the descriptions of reads in a FASTQ file."""
    matches_per_read = []
    
    # Handle input as stdin when fastq_file is '-'
    if fastq_file == '-':
        print("Reading from standard input...")
        fastq_handle = sys.stdin
    else:
        fastq_handle = open(fastq_file, "r")
    
    for record in SeqIO.parse(fastq_handle, "fastq"):
        # Extract the number of matches from the description
        desc = record.description
        match_info = re.search(r'matches=(\d+)', desc)
        if match_info:
            matches_per_read.append(int(match_info.group(1)))
    
    if fastq_file != '-':
        fastq_handle.close()
    
    return matches_per_read

def store_results(classified_reads, output_file):
    """Store the classified reads in a FASTQ file."""
    with open(output_file, 'w') as out_file:
        SeqIO.write(classified_reads, out_file, "fastq")

def count_reads(fastq_file):
    """Count the number of reads in a FASTQ file."""
    # Handle input as stdin when fastq_file is '-'
    if fastq_file == '-':
        # For stdin, we can't count without consuming, so return -1
        print("Reading from stdin, cannot count reads in advance")
        return -1
    else:
        return sum(1 for _ in SeqIO.parse(fastq_file, "fastq"))
    


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description='Classify FASTQ reads based on reference sequences')
    parser.add_argument('--ref', '-r', required=True, help='Reference file with SNPs and indels')
    parser.add_argument('--fastq', '-f', required=True, help='Input FASTQ file to process or - for stdin')
    parser.add_argument('--output', '-o', required=True, help='Output file for classified reads')
    parser.add_argument('--plot', '-p', default=None, help='Output file for histogram plot (optional)')
    parser.add_argument('--threshold', '-t', type=int, default=15, help='Match threshold (default: 15)')
    parser.add_argument('--processes', '-n', type=int, default=None, help='Number of processes to use (default: all available CPU cores)')
    parser.add_argument('--histogram_only', action='store_true', help='Only generate histogram from already processed reads')
    parser.add_argument('--stats', '-s', default=None, help='Output file for match statistics (one count per line)')
    parser.add_argument('--reverse_complement', '-rc', action='store_true', 
                      help='Include reverse complements of reference k-mers (default: True)')
    
    args = parser.parse_args()
    
    # Use provided arguments
    reference_file = args.ref
    fastq_file = args.fastq
    output_file = args.output
    plot_histo = args.plot is not None


    # Wenn nur Histogramm erstellt werden soll
    if args.histogram_only:
        print("Running in histogram-only mode...")
        print(f"Extracting match counts from processed reads...")
        matches_per_read = extract_match_counts(fastq_file)
        if matches_per_read:
            plot_histogram(matches_per_read, args.threshold, args.plot)
        else:
            print("No match information found in the reads. Make sure the reads were processed with this script.")
    else:
        # Normal operation - count, classify, and store reads
        print(f"Counting reads in input file...")
        initial_read_count = count_reads(fastq_file)
        if initial_read_count > 0:
            print(f"Initial read count: {initial_read_count}")

        print(f"Generating reference hashes...")
        reference_hashes, removed_count = generate_reference_hashes(reference_file, include_reverse_complement=not args.reverse_complement)
        print(f"Generated {len(reference_hashes)} reference hashes ({removed_count} sequences skipped)")

        print(f"Classifying reads with threshold {args.threshold}...")
        classified_reads = classify_reads(fastq_file, reference_hashes, threshold=args.threshold, 
                                        plot_histo=plot_histo, num_processes=args.processes, stats_file=args.stats)
        
        print(f"Storing classified reads...")
        store_results(classified_reads, output_file)
        
     
        # Count the number of reads in the output FASTQ file
        final_read_count = len(classified_reads)
        print(f"Final read count: {final_read_count}")
        if initial_read_count > 0:
            print(f"Classification complete. {final_read_count} of {initial_read_count} reads classified.")
        else:
            print(f"Classification complete. {final_read_count} reads classified.")
