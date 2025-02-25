## This script will take the concatenated library of all TE detected from the pangenome and curatedlib TE and will process it into a non-redundant quality panTElib

import argparse
import glob, os
import subprocess
from Bio import SeqIO

parser=argparse.ArgumentParser(description="Generate panTE library")
parser.add_argument('--all_TElib', type=str, help='Path to concatenated TElib')
parser.add_argument('--curatedlib', type=str, help='Path to curated library')
parser.add_argument('--threads', type=str, help='Number of threads to be used')
args=parser.parse_args()

# Input Files
all_TElib = args.all_TElib
curatedlib = args.curatedlib
threads = args.threads

# Script dir
snake_dir = os.getcwd()

# Step 1: Run TEsorter and postprocess
def generate_fasta_length_table(input_file):
    """
    Generates a table with sequence IDs and lengths from a FASTA file.
    :param input_file: Input FASTA file.
    :return: Dictionary with sequence IDs as keys and lengths as values.
    """
    fasta_table = {}
    sequence_id = None
    sequence_length = 0

    with open(input_file, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                if sequence_id:
                    fasta_table[sequence_id] = sequence_length
                sequence_id = line[1:].strip().split()[0]  # Extract ID from header
                sequence_length = 0
            else:
                sequence_length += len(line.strip())

        if sequence_id:
            fasta_table[sequence_id] = sequence_length

    return fasta_table

# Here we run TEsorter to find the TE families that have enough sequence information to be correctly classified
def run_tesorter_and_postprocess(db, threads, input_file):
    """
    Runs TEsorter with the specified parameters and processes the output.
    :param db: The database to use with TEsorter (e.g., 'rexdb-plant').
    :param threads: Number of threads to use for TEsorter.
    :param input_file: Input file path for TEsorter.
    """
    # Run TEsorter
    tesorter_command = f"TEsorter -db {db} -p {threads} -dp2 {input_file}"
    subprocess.run(tesorter_command, shell=True, check=True)

    # Process the output
    cls_file = f"{input_file}.{db}.cls.tsv"
    processed_list_file = f"{input_file}.classified.list"

    with open(cls_file, 'r') as infile, open(processed_list_file, 'w') as outfile:
        next(infile)  # Skip the header line
        for line in infile:
            processed_line = line.split("\t")[0]
            outfile.write(f"{processed_line}\n")

    # Generate FASTA table
    fasta_table = generate_fasta_length_table(input_file)

    # Generate classified and unclassified outputs
    classified_list = set()
    with open(processed_list_file, 'r') as classified:
        for line in classified:
            classified_list.add(line.strip().split('\t')[0])

    final_output = f"{input_file}.table.classunclass"
    with open(final_output, 'w') as output:
        output.write("genome\tte_name\torder\tsuperfamily\tlength\tcl_uncl\n")  # Write header
        for seq_id, seq_length in fasta_table.items():
            parts = seq_id.split("_TE_")
            genome = parts[0]
            te_name_order_superfamily = parts[1].split("#") if len(parts) > 1 else ["unknown", "unknown"]
            te_name = "TE_"+te_name_order_superfamily[0]
            order_superfamily = te_name_order_superfamily[1].split("/") if len(te_name_order_superfamily) > 1 else ["unknown", "unknown"]
            order = order_superfamily[0]
            superfamily = order_superfamily[1] if len(order_superfamily) > 1 else "unknown"
            classification = "classified" if seq_id in classified_list else "unclassified"

            output.write(f"{genome}\t{te_name}\t{order}\t{superfamily}\t{seq_length}\t{classification}\n")


# Step 2:
# extract TEs that pass the filter and sort by length in descending order

def extract_filtered_sequences_from_allTElib(input_file, fasta_file, output_file):
    # Step 1: Read the input file and extract TE names, order, and superfamily
    t_data = []
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        for line in lines[1:]:  # Skipping the header
            line = line.replace('"','')
            cols = line.strip().split(' ')
            genome = cols[1]
            te_name = cols[2]
            order = cols[3]
            superfamily = cols[4]
            
            # Step 2: Create the formatted ID for the FASTA file
            fasta_id = f"{genome}_{te_name}#{order}/{superfamily}"
            
            t_data.append(fasta_id)

    # Step 3: Extract sequences from FASTA file that match the IDs
    matched_sequences = []
    
    with open(fasta_file, 'r') as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id in t_data:
                matched_sequences.append(record)
    
    # Step 4: Sort the sequences by length in descending order
    matched_sequences.sort(key=lambda seq: len(seq), reverse=True)
    
    # Step 5: Write the sorted sequences to the output FASTA file
    with open(output_file, 'w') as output_fasta:
        SeqIO.write(matched_sequences, output_fasta, "fasta")
    
    

def concatenate_fasta_files(fasta_file1, fasta_file2, output_file):
    # Read sequences from both FASTA files
    sequences = []
    
    # Read first FASTA file
    with open(fasta_file1, 'r') as file1:
        sequences.extend(SeqIO.parse(file1, "fasta"))
    
    # Read second FASTA file
    with open(fasta_file2, 'r') as file2:
        sequences.extend(SeqIO.parse(file2, "fasta"))
    
    # Write concatenated sequences to the output file
    with open(output_file, 'w') as output_fasta:
        SeqIO.write(sequences, output_fasta, "fasta")



# Main function
def main():

    # go to allTElib directory (necessary for TEsorter to work fine)
    dir=os.path.dirname(all_TElib)
    if dir!='':
        os.chdir(dir)
    
    # Step 1: TEsorter get classification and postprocess
    run_tesorter_and_postprocess('rexdb-plant', threads, all_TElib)

    # Step 2: Filter based on lengths
    R_scipt_run=f"Rscript {snake_dir}/scripts_for_snakepipe/lenTEandprotfilter.R -i {all_TElib}.table.classunclass" 
    subprocess.run(R_scipt_run, shell=True, check=True)

    # Step 3: Extracted filter pass from allTElib
    extraction_input=f"{all_TElib}.table.classunclass.filtered_TEs.txt"
    extraction_output=f"{all_TElib}.lenprotfiltr.fa"
    extract_filtered_sequences_from_allTElib(extraction_input, all_TElib, extraction_output)

    # Step 4: Concatenate with curatedlib and run vsearch
    curatedlib_name=os.path.splitext(os.path.basename(curatedlib))[0]
    concat_output=f"{all_TElib}.{curatedlib_name}.fa"
    concatenate_fasta_files(curatedlib, extraction_output, concat_output)

    vsearch_output="panTElib.vsearch.centroids.fa"
    blast_output="panTElib.vsearch.results.blast6"
    vsearch_run=f"vsearch --cluster_smallmem {concat_output} --usersort --id 0.8 --query_cov 0.8 --target_cov 0.8 --threads {threads} --centroids {vsearch_output} --blast6out {blast_output}" 
    subprocess.run(vsearch_run, shell=True, check=True)

    # remove temp
    for f in glob.glob("*rexdb-plant*"):
        os.remove(f)

# Run the Script
if __name__ == "__main__":
    main()