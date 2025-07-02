import os
import sys
import csv
import datetime
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter
import mappy as mp  # minimap2 for alignment
from Bio.Align.Applications import MafftCommandline
import tempfile
import subprocess

# Check sys.arg
if len(sys.argv) < 4:
    print("Error: Please provide the reference directory, minimum read length, and maximum read length.")
    exit()

# --- INPUTS ---
BARCODE_FASTA = f"{sys.argv[1]}/barcodes.fasta"
REFERENCE_FASTA = f"{sys.argv[1]}/reference.fasta"
FASTQ_FILE = f"{sys.argv[1]}/input_reads.fastq"
LENGTH_CUTOFF = (int(sys.argv[2]), int(sys.argv[3])) # 2000 and 3000 for MPH

# --- OUTPUT FOLDER ---
today = datetime.datetime.now().strftime("%y%m%d")
OUTPUT_DIR = f"{sys.argv[1]}/extracted_fastq_{today}"
SUMMARY_CSV = f"{sys.argv[1]}/demultiplex_summary.csv"

# --- PARAMETERS ---
MIN_ALIGNMENT_IDENTITY = 0.8  # Minimum identity for a read to be considered a hit

# --- 1. Parse barcodes ---
def parse_barcodes(barcode_fasta):
    fwd_barcodes = {}
    rev_barcodes = {}
    for record in SeqIO.parse(barcode_fasta, "fasta"):
        if record.id.lower().startswith("fwd"):
            fwd_barcodes[record.id] = str(record.seq).upper()
        elif record.id.lower().startswith("rev"):
            rev_barcodes[record.id] = str(record.seq).upper()
    return fwd_barcodes, rev_barcodes

# --- 2. Plate layout mapping ---
def well_name(row, col):
    return f"{chr(65+row)}{col+1}"

def get_well(fwd_idx, rev_idx):
    # fwd_idx: 0-11, rev_idx: 0-7
    return well_name(rev_idx, fwd_idx)

# --- 3. Demultiplex reads ---
def find_barcodes_in_read(seq, fwd_barcodes, rev_barcodes):
    # Returns (fwd_id, rev_id, start, end, orientation)
    seq_upper = seq.upper()
    for fwd_id, fwd_seq in fwd_barcodes.items():
        fwd_pos = seq_upper.find(fwd_seq)
        if fwd_pos != -1:
            for rev_id, rev_seq in rev_barcodes.items():
                rev_rc = str(Seq(rev_seq).reverse_complement())
                rev_rc_pos = seq_upper.find(rev_rc)
                if rev_rc_pos != -1 and rev_rc_pos > fwd_pos + len(fwd_seq):
                    # Forward orientation
                    start = fwd_pos + len(fwd_seq)
                    end = rev_rc_pos
                    return (fwd_id, rev_id, start, end, "fwd")
    # Try reverse orientation
    for rev_id, rev_seq in rev_barcodes.items():
        rev_pos = seq_upper.find(rev_seq)
        if rev_pos != -1:
            for fwd_id, fwd_seq in fwd_barcodes.items():
                fwd_rc = str(Seq(fwd_seq).reverse_complement())
                fwd_rc_pos = seq_upper.find(fwd_rc)
                if fwd_rc_pos != -1 and fwd_rc_pos > rev_pos + len(rev_seq):
                    # Reverse orientation
                    start = rev_pos + len(rev_seq)
                    end = fwd_rc_pos
                    return (fwd_id, rev_id, start, end, "rev")
    return None

# --- 4. Align to reference ---
def align_to_reference(seq, aligner):
    # Returns (is_significant, ref_start, ref_end, identity)
    for hit in aligner.map(seq):
        # Best hit with high quality
        if hit.is_primary and hit.mapq > 0:
            identity = 1 - (hit.NM / hit.blen)
            if identity >= MIN_ALIGNMENT_IDENTITY:
                return True, hit.r_st, hit.r_en, identity
    return False, None, None, None

# --- 5. Robust consensus and mutation calling ---
def robust_consensus_and_mutations(filtered_reads, reference_seq, output_dir, well_name, threshold=0.5):
    """
    filtered_reads: list of str (filtered reads for a well)
    reference_seq: str (reference sequence)
    output_dir: str (directory to save alignment file)
    well_name: str (well identifier for file naming)
    threshold: float (fraction of sequences with consensus for base calling)
    Returns: consensus_seq (str), mutations (list of str)
    """
    if len(filtered_reads) == 0:
        return None, []
    
    # Add reference as the first sequence
    all_seqs = [reference_seq] + filtered_reads
    
    # Create alignment file path in the well's directory
    alignment_file = os.path.join(output_dir, f"{well_name}_alignment.fasta")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        reads_fasta = os.path.join(tmpdir, "reads.fasta")
        aligned_fasta = os.path.join(tmpdir, "aligned.fasta")
        SeqIO.write(
            [SeqRecord(Seq(s), id=f"seq{i}") for i, s in enumerate(all_seqs)],
            reads_fasta,
            "fasta"
        )
        mafft_cline = MafftCommandline(input=reads_fasta)
        with open(aligned_fasta, "w") as handle:
            subprocess.run(str(mafft_cline), shell=True, stdout=handle, check=True)
        
        # Copy alignment file to the well's directory
        with open(aligned_fasta, "r") as src, open(alignment_file, "w") as dst:
            dst.write(src.read())
        
        alignment = AlignIO.read(aligned_fasta, "fasta")
        ref_aligned = str(alignment[0].seq)
        read_aligned = [str(rec.seq) for rec in alignment[1:]]

        consensus_seq = ""
        mutations = []
        for i in range(len(ref_aligned)):
            ref_base = ref_aligned[i]
            read_bases = [r[i] for r in read_aligned]
            base_counts = Counter(read_bases)
            # Do not exclude gaps from consensus calculation
            base_counts_no_N = {b: c for b, c in base_counts.items() if b != "N"}
            if not base_counts_no_N:
                consensus_seq += "N"
                continue
            most_common_base, most_common_count = Counter(base_counts_no_N).most_common(1)[0]
            total_reads = sum(base_counts_no_N.values())
            # Majority of reads differ from reference and agree (including "-")
            if most_common_base != ref_base and most_common_count / total_reads >= threshold:
                consensus_seq += most_common_base.upper()
                mutations.append(f"{ref_base.upper()}{i+1}{most_common_base.upper()}")
            # Majority of reads agree with reference
            elif ref_base in base_counts_no_N and base_counts_no_N[ref_base] / total_reads >= threshold:
                consensus_seq += ref_base.upper()
            else:
                consensus_seq += "N"
        if not mutations:
            mutations.append("WT")
        return consensus_seq, mutations

# --- 6. Main pipeline ---
def main():
    # 1. Parse barcodes
    fwd_barcodes, rev_barcodes = parse_barcodes(BARCODE_FASTA)
    fwd_list = sorted(fwd_barcodes.keys())
    rev_list = sorted(rev_barcodes.keys())

    # 2. Load reference
    ref_record = next(SeqIO.parse(REFERENCE_FASTA, "fasta"))
    ref_seq = str(ref_record.seq).upper()
    aligner = mp.Aligner(REFERENCE_FASTA, preset="map-ont")

    # 3. Prepare output folders
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    well_fastq_handles = {}
    for r in range(8):
        for c in range(12):
            well = well_name(r, c)
            well_dir = os.path.join(OUTPUT_DIR, well)
            os.makedirs(well_dir, exist_ok=True)
            fq_path = os.path.join(well_dir, f"{well}.fastq")
            well_fastq_handles[well] = open(fq_path, "w")

    # 4. Demultiplex and align
    well_reads = defaultdict(list)
    for record in SeqIO.parse(FASTQ_FILE, "fastq"):
        if len(str(record.seq)) < LENGTH_CUTOFF[0] or len(str(record.seq)) > LENGTH_CUTOFF[1]:
            continue
        result = find_barcodes_in_read(str(record.seq), fwd_barcodes, rev_barcodes)
        if result:
            fwd_id, rev_id, start, end, orientation = result
            sub_seq = str(record.seq)[start:end]
            if orientation == "rev":
                sub_seq = str(Seq(sub_seq).reverse_complement())
            is_sig, _, _, identity = align_to_reference(sub_seq, aligner)
            if is_sig:
                fwd_idx = fwd_list.index(fwd_id)
                rev_idx = rev_list.index(rev_id)
                well = get_well(fwd_idx, rev_idx)
                SeqIO.write(record, well_fastq_handles[well], "fastq")
                well_reads[well].append(sub_seq)

    # Close all FASTQ handles
    for handle in well_fastq_handles.values():
        handle.close()

    # 5. Consensus and mutation calling
    with open(SUMMARY_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Well", "Consensus", "Mutations", "ReadCount"])
        for well in sorted(well_reads.keys()):
            reads = well_reads[well]
            if not reads:
                continue
            well_output_dir = os.path.join(OUTPUT_DIR, well)
            consensus_seq, mutations = robust_consensus_and_mutations(reads, ref_seq, well_output_dir, well)
            writer.writerow([well, consensus_seq, ";".join(mutations), len(reads)])


if __name__ == "__main__":
    main()