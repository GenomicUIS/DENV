from Bio import SeqIO
# set max_ambigous to your specificities 
def is_valid_sequence(seq, max_ambiguous=0): 
    """
    Check if a sequence is valid based on the fraction of ambiguous bases.
    
    Parameters:
        seq (str): The sequence to check.
        max_ambiguous (float): Maximum allowed fraction of ambiguous bases (default: 5%).
    
    Returns:
        bool: True if the sequence is valid, False otherwise.
    """
    # Count ambiguous bases (non-ACGT characters)
    ambiguous_count = sum(1 for base in seq if base.upper() not in "ACGT")
    
    # Calculate the fraction of ambiguous bases
    fraction_ambiguous = ambiguous_count / len(seq)
    
    # Return True if the fraction is within the allowed limit
    return fraction_ambiguous <= max_ambiguous

def filter_fasta(input_file, output_file, max_ambiguous=0):
    """
    Filter sequences in a FASTA file based on the presence of ambiguous bases.
    
    Parameters: 
        input_file (str): "Replace with your input file path"
        output_file (str): "Replace with your input file path"
        max_ambiguous (float): Maximum allowed fraction of ambiguous bases (default: 5%).
    """
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            if is_valid_sequence(record.seq, max_ambiguous):
                SeqIO.write(record, out_handle, "fasta")

# File paths
input_file = "sequence.FASTA"  # Replace with your input file path
output_file = "output.FASTA"  # Replace with your desired output file path

# Filter the FASTA file
filter_fasta(input_file, output_file)

print(f"Filtered sequences saved to {output_file}")