import sys
import os
import csv
import datetime
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
import tempfile
import subprocess
import pandas as pd

# Check sys.arg
if len(sys.argv) < 2:
    print("Error: Please provide the reference directory.")
    exit()

# --- INPUTS ---
DEMULTIPLEX_CSV = f"{sys.argv[1]}/demultiplex_summary.csv"
REFERENCE_ORF_FASTA = f"{sys.argv[1]}/reference_orf.fasta"

# --- OUTPUT ---
today = datetime.datetime.now().strftime("%y%m%d")
AA_MUTATIONS_CSV = f"{sys.argv[1]}/amino_acid_mutations_{today}.csv"

def align_dna_to_orf(consensus_dna_seqs, ref_orf_seq, well_names):
    """
    Align consensus DNA sequences to reference ORF using MAFFT
    Returns aligned sequences with ORF region extracted
    """
    if not consensus_dna_seqs:
        return None, {}
    
    # Create sequence records
    all_seqs = [SeqRecord(Seq(ref_orf_seq), id="reference_orf")]
    for i, (dna_seq, well) in enumerate(zip(consensus_dna_seqs, well_names)):
        if dna_seq:  # Only add non-empty sequences
            all_seqs.append(SeqRecord(Seq(dna_seq), id=well))
    
    if len(all_seqs) <= 1:
        return None, {}
    
    with tempfile.TemporaryDirectory() as tmpdir:
        dna_fasta = os.path.join(tmpdir, "dna_seqs.fasta")
        aligned_fasta = os.path.join(tmpdir, "aligned_dna.fasta")
        
        # Write sequences to file
        SeqIO.write(all_seqs, dna_fasta, "fasta")
        
        # Run MAFFT alignment
        mafft_cline = MafftCommandline(input=dna_fasta)
        with open(aligned_fasta, "w") as handle:
            subprocess.run(str(mafft_cline), shell=True, stdout=handle, check=True)
        
        # Read alignment
        alignment = AlignIO.read(aligned_fasta, "fasta")
        ref_orf_aligned = str(alignment[0].seq)
        consensus_aligned = {}
        
        for record in alignment[1:]:
            consensus_aligned[record.id] = str(record.seq)
        
        return ref_orf_aligned, consensus_aligned

def extract_orf_region_and_translate(ref_orf_aligned, consensus_aligned, well_name):
    """
    Extract ORF region from aligned consensus and translate to amino acids
    Uses alignment to reference ORF to determine the correct region
    Preserves alignment information for proper deletion detection
    """
    if well_name not in consensus_aligned:
        return "", "No_consensus"
    
    consensus_seq = consensus_aligned[well_name]
    
    if len(ref_orf_aligned) != len(consensus_seq):
        return "", "Alignment_error"
    
    # Extract ORF region by removing gaps that correspond to gaps in reference
    orf_consensus = ""
    for i in range(len(ref_orf_aligned)):
        ref_base = ref_orf_aligned[i]
        cons_base = consensus_seq[i]
        
        # Only keep positions where reference has a base (not a gap)
        if ref_base != '-':
            if cons_base == '-':
                # Gap in consensus where reference has base - this is a deletion
                # Use a special marker that will translate to a gap in protein
                orf_consensus += '---'  # Triple dash to maintain reading frame
            else:
                orf_consensus += cons_base
    
    # For translation, we need to handle the deletion markers
    # Split into codons and translate each, preserving gaps
    protein_with_gaps = ""
    for i in range(0, len(orf_consensus), 3):
        codon = orf_consensus[i:i+3]
        if len(codon) == 3:
            if codon == '---':
                protein_with_gaps += '-'  # Deletion in protein
            elif 'N' in codon or '-' in codon:
                protein_with_gaps += 'X'  # Ambiguous amino acid
            else:
                try:
                    aa = str(Seq(codon).translate())
                    protein_with_gaps += aa
                except:
                    protein_with_gaps += 'X'
        else:
            # Incomplete codon at end
            break
    
    return protein_with_gaps, orf_consensus

def translate_dna(dna_seq):
    """
    Simple DNA to protein translation
    """
    if len(dna_seq) < 3:
        return ""
    
    # Ensure sequence length is multiple of 3
    while len(dna_seq) % 3 != 0:
        dna_seq = dna_seq[:-1]
    
    try:
        protein = str(Seq(dna_seq).translate())
        return protein
    except Exception as e:
        print(f"Translation error: {e}")
        return ""

def align_protein_sequences(ref_protein, consensus_proteins, well_names):
    """
    Align reference protein with consensus proteins using MAFFT
    Returns aligned sequences
    """
    if not consensus_proteins:
        return None, []
    
    # Create sequence records
    all_seqs = [SeqRecord(Seq(ref_protein), id="reference")]
    for i, (protein, well) in enumerate(zip(consensus_proteins, well_names)):
        if protein:  # Only add non-empty proteins
            all_seqs.append(SeqRecord(Seq(protein), id=well))
    
    if len(all_seqs) <= 1:
        return None, []
    
    with tempfile.TemporaryDirectory() as tmpdir:
        protein_fasta = os.path.join(tmpdir, "proteins.fasta")
        aligned_fasta = os.path.join(tmpdir, "aligned_proteins.fasta")
        
        # Write sequences to file
        SeqIO.write(all_seqs, protein_fasta, "fasta")
        
        # Run MAFFT alignment
        mafft_cline = MafftCommandline(input=protein_fasta)
        with open(aligned_fasta, "w") as handle:
            subprocess.run(str(mafft_cline), shell=True, stdout=handle, check=True)
        
        # Read alignment
        alignment = AlignIO.read(aligned_fasta, "fasta")
        ref_aligned = str(alignment[0].seq)
        consensus_aligned = {}
        
        for i, record in enumerate(alignment[1:]):
            consensus_aligned[record.id] = str(record.seq)
        
        return ref_aligned, consensus_aligned

def find_amino_acid_mutations(ref_aligned, consensus_aligned, well_name):
    """
    Compare aligned reference and consensus proteins to find mutations
    Returns list of mutation strings (e.g., "A123V", "S193del")
    """
    mutations = []
    
    if well_name not in consensus_aligned:
        return ["No_consensus"]
    
    consensus_seq = consensus_aligned[well_name]
    
    if len(ref_aligned) != len(consensus_seq):
        return ["Alignment_error"]
    
    # Track position in reference sequence (ignoring gaps)
    ref_pos = 0
    
    for i in range(len(ref_aligned)):
        ref_aa = ref_aligned[i]
        cons_aa = consensus_seq[i]
        
        # Only increment position for non-gap positions in reference
        if ref_aa != '-':
            ref_pos += 1
            
            # Check for deletions (gap in consensus where reference has amino acid)
            if cons_aa == '-':
                mutations.append(f"{ref_aa}{ref_pos}del")
            # Check for substitutions (ignoring ambiguous amino acids)
            elif cons_aa != 'X' and ref_aa != cons_aa:
                mutations.append(f"{ref_aa}{ref_pos}{cons_aa}")
        # Check for insertions (amino acid in consensus where reference has gap)
        elif ref_aa == '-' and cons_aa != '-' and cons_aa != 'X':
            # For insertions, we use the previous reference position
            mutations.append(f"ins{ref_pos}{cons_aa}")
    
    return mutations if mutations else ["WT"]

def main():
    print("Starting amino acid analysis...")
    
    # 1. Load reference ORF and translate
    print("Loading reference ORF...")
    try:
        ref_orf_record = next(SeqIO.parse(REFERENCE_ORF_FASTA, "fasta"))
        ref_orf_seq = str(ref_orf_record.seq).upper()
        ref_protein = translate_dna(ref_orf_seq)
        print(f"Reference ORF length: {len(ref_orf_seq)} bp")
        print(f"Reference protein length: {len(ref_protein)} aa")
    except Exception as e:
        print(f"Error loading reference ORF: {e}")
        return
    
    # 2. Load demultiplex summary
    print("Loading demultiplex summary...")
    try:
        df = pd.read_csv(DEMULTIPLEX_CSV)
        print(f"Loaded {len(df)} wells from demultiplex summary")
    except Exception as e:
        print(f"Error loading demultiplex CSV: {e}")
        return
    
    # 3. Align consensus DNA sequences to reference ORF
    print("Aligning consensus DNA to reference ORF...")
    results = []
    consensus_dna_seqs = []
    well_names = []
    
    # Collect consensus sequences for alignment
    for _, row in df.iterrows():
        well = row['Well']
        consensus_dna = row['Consensus']
        read_count = row['ReadCount']
        
        if pd.isna(consensus_dna) or consensus_dna == "":
            results.append({
                'Well': well,
                'ReadCount': read_count,
                'ORF_DNA': "",
                'ConsensusProtein': "",
                'AAMutations': "No_consensus",
                'NumAAMutations': 0
            })
            continue
        
        # Keep original consensus (with gaps) for alignment
        consensus_dna_seqs.append(consensus_dna)
        well_names.append(well)
        
        # Initialize result entry
        results.append({
            'Well': well,
            'ReadCount': read_count,
            'ORF_DNA': "",
            'ConsensusProtein': "",
            'AAMutations': "",
            'NumAAMutations': 0
        })
    
    # 4. Align DNA sequences to ORF and extract ORF regions
    print("Extracting ORF regions and translating...")
    if consensus_dna_seqs:
        ref_orf_aligned, consensus_dna_aligned = align_dna_to_orf(consensus_dna_seqs, ref_orf_seq, well_names)
        
        if ref_orf_aligned and consensus_dna_aligned:
            # Extract ORF regions and translate
            consensus_proteins = []
            orf_well_names = []
            
            for result in results:
                well = result['Well']
                if well in consensus_dna_aligned:
                    consensus_protein, orf_dna = extract_orf_region_and_translate(ref_orf_aligned, consensus_dna_aligned, well)
                    result['ORF_DNA'] = orf_dna
                    result['ConsensusProtein'] = consensus_protein
                    
                    if consensus_protein and not consensus_protein.startswith("Translation_error") and consensus_protein != "Too_short":
                        consensus_proteins.append(consensus_protein)
                        orf_well_names.append(well)
                    else:
                        result['AAMutations'] = consensus_protein if consensus_protein else "Translation_failed"
                else:
                    result['AAMutations'] = "DNA_alignment_failed"
    
    # 5. Align proteins and find mutations
    print("Aligning proteins and finding mutations...")
    if consensus_proteins:
        ref_aligned, protein_aligned = align_protein_sequences(ref_protein, consensus_proteins, orf_well_names)
        
        if ref_aligned and protein_aligned:
            # Update results with mutation information
            for result in results:
                well = result['Well']
                if well in protein_aligned and result['AAMutations'] == "":
                    mutations = find_amino_acid_mutations(ref_aligned, protein_aligned, well)
                    result['AAMutations'] = ";".join(mutations)
                    result['NumAAMutations'] = len(mutations) if mutations != ["WT"] else 0
    
    # 6. Export results
    print("Exporting results...")
    with open(AA_MUTATIONS_CSV, 'w', newline='') as csvfile:
        fieldnames = ['Well', 'ReadCount', 'ORF_DNA', 'ConsensusProtein', 'AAMutations', 'NumAAMutations']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Results exported to {AA_MUTATIONS_CSV}")
    
    # Print summary statistics
    total_wells = len(results)
    wells_with_reads = sum(1 for r in results if r['ReadCount'] > 0)
    wells_with_protein = sum(1 for r in results if r['ConsensusProtein'] != "")
    wells_with_mutations = sum(1 for r in results if r['NumAAMutations'] > 0)
    
    print(f"\nSummary:")
    print(f"Total wells: {total_wells}")
    print(f"Wells with reads: {wells_with_reads}")
    print(f"Wells with translatable protein: {wells_with_protein}")
    print(f"Wells with amino acid mutations: {wells_with_mutations}")

if __name__ == "__main__":
    main()
