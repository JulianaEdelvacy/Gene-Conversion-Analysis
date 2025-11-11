from Bio.SeqIO import parse, write
from Bio import pairwise2
from statistics import mean
import sys
from io import TextIOWrapper


def check_read_quality(assembled_read, R1_ids_list, R2_ids_list, R1_seqs, R2_seqs, output_file: TextIOWrapper):
    """
    Checks the quality of the assembled read's regions (non-overlapping R1, 
    non-overlapping R2, and the overlapping middle region) based on phred scores 
    and writes the read to the output file if quality passes thresholds.
    
    Args:
        assembled_read (SeqRecord): The assembled sequence record (with quality scores).
        R1_ids_list (list): List of sequence IDs from the R1 file.
        R2_ids_list (list): List of sequence IDs from the R2 file.
        R1_seqs (list): List of sequence objects from the R1 file.
        R2_seqs (list): List of sequence objects from the R2 file.
        output_file (TextIOWrapper): File handle for writing selected sequences.
    """
    
    read_id = assembled_read.id
    assembled_seq = assembled_read.seq
    seq_len = len(assembled_seq)

    # --- R1 ALIGNMENT (Assuming R1 starts at the beginning of the assembly) ---
    try:
        # Get the index of the assembled sequence ID from R1's list of IDs
        index_temp = R1_ids_list.index(read_id)
    except ValueError:
        return None # Assembled read ID not found in R1 list

    R1_seq = R1_seqs[index_temp]
    R1_len = len(R1_seq)

    # Target sequence for R1 alignment (from assembly start)
    target_seq_R1 = assembled_seq[0:R1_len]
    R1_alignment_test = pairwise2.align.localxx(R1_seq, target_seq_R1)
    
    if len(R1_alignment_test) == 0:
         return None
         
    R1_local_assembly_alignment = R1_alignment_test[0]


    # --- R2 ALIGNMENT (Starts near the original snippet logic) ---
    
    try:
        # Get the index of the assembled sequence ID from R2's list of IDs
        index_temp = R2_ids_list.index(read_id) 
    except ValueError:
        return None # Assembled read ID not found in R2 list

    R2_seq = R2_seqs[index_temp]
    
    # Index of probable start of R2 in the assembly
    R2_start_index = seq_len - len(R2_seq) 
    
    # Target sequence for R2 alignment (from assembly end)
    target_seq_R2 = assembled_seq[R2_start_index:seq_len] 
    
    # Align R2's reverse complement to the target sequence in the assembly
    R2_alignment_test = pairwise2.align.localxx(R2_seq.reverse_complement(), target_seq_R2)   
    
    if len(R2_alignment_test) == 0:
        return None
        
    R2_local_assembly_alignment = R2_alignment_test[0]
    R2_len = len(R2_seq)

    # Define the intersection/overlap region of R1 and R2 reads relative to the assembly
    # This range is from where R2 begins (R2_start_index) up to where R1 ends (R1_len)
    overlap_range = range(R2_start_index, R1_len)
    
    '''
    Alignment Conditions Check:
    
    The logic ensures that the gap-free alignment region intersects with the 
    expected read overlap region.
    
    1. R1's gap-free region starts at 0. Check if the END of this region 
       falls within the overlap_range.
    2. R2's gap-free region ends at the assembly end. Check if the START 
       of this region falls within the overlap_range.
    '''
    
    # R1 region that aligns with assembly without gaps (from the start of the alignment string)
    R1_no_gap_region = R1_local_assembly_alignment[0].split("-")[0] 
    
    # Check if the last index of the R1 non-gap region is inside the overlap range
    if (len(R1_no_gap_region) - 1) not in overlap_range:
        return None
        
    # R2 region that aligns with assembly without gaps (from the end of the alignment string)
    R2_no_gap_region = R2_local_assembly_alignment[0].split("-")[-1] 
    
    # Calculate the starting index of the R2 non-gap region in the assembly
    R2_no_gap_start_index = seq_len - len(R2_no_gap_region) 
    
    # Check if the first index of the R2 non-gap region is inside the overlap range
    if (R2_no_gap_start_index - 1) not in overlap_range:
        return None

    # --- QUALITY SCORE FILTERING ---
    
    # 1. Initial Region (R1 non-overlapping part): 
    #    Starts at index 0 and ends where R2 begins (R2_start_index)
    initial_region_qual = assembled_read.letter_annotations["phred_quality"][0:R2_start_index] 
    
    # 2. Final Region (R2 non-overlapping part): 
    #    Starts where R1 ends (R1_len) and ends at the sequence end (seq_len)
    final_region_qual = assembled_read.letter_annotations["phred_quality"][R1_len:seq_len] 
    
    # 3. Middle Region (R1/R2 Overlapping part): 
    #    Starts where R2 begins (R2_start_index) and ends where R1 ends (R1_len)
    middle_region_qual = assembled_read.letter_annotations["phred_quality"][R2_start_index:R1_len] 
    
    try:
        # Quality Threshold Conditions:
        # Minimum quality in non-overlapping regions > 20
        # Mean quality in the overlapping region > 30
        if (min(initial_region_qual) > 20) and \
           (min(final_region_qual) > 20) and \
           (mean(middle_region_qual) > 30):   
           
            # Write the assembled read to the output file if it passes the quality check
            write(assembled_read, output_file, "fastq")
            
    except ValueError:
        # Catches error if a region is empty (e.g., min() of empty sequence)
        return None
        
    except Exception as e:
        # Catch other potential errors
        print(f"Error processing read {read_id}: {e}", file=sys.stderr)
        return None
        
    return None

# --- MAIN SCRIPT EXECUTION ---

# List of assembled sequences (SeqRecord objects from the assembled FASTQ file)
assembled_seqs = list(parse("assemble-pass_test.fastq", "fastq")) 

# List of R1 sequences
R1_reads = list(parse("R1_test.fastq", "fastq")) 
# List of R2 sequences
R2_reads = list(parse("R2_test.fastq", "fastq")) 

# Extracting ID and Sequence lists for quick lookup (as done in original script)
R1_ids_list = [read.id for read in R1_reads]
R1_seqs = [read.seq for read in R1_reads]

R2_ids_list = [read.id for read in R2_reads]
R2_seqs = [read.seq for read in R2_reads]

# Open the output file in append mode
selected_seqs_file = open("selected_seqs.fastq", "a+")

print(f"Starting quality check on {len(assembled_seqs)} assembled reads...")

for i in range(len(assembled_seqs)):
    
    # Re-open the file periodically to flush the buffer (every 1000 reads)
    if (i % 1000 == 0) and (i != 0):
        print(f"Processing read {i}")
        selected_seqs_file.close()
        selected_seqs_file = open("selected_seqs.fastq", "a+")
        
    # Check the quality of the current assembled read
    check_read_quality(assembled_seqs[i], R1_ids_list, R2_ids_list, R1_seqs, R2_seqs, selected_seqs_file)

# Close the final file handle
selected_seqs_file.close()

print("Quality check complete. Selected sequences saved to selected_seqs.fastq.")
