# usage: python3 add_fullsizeLTRRTs.py EDTA_TElib_fa EDTA_intact_fa work_dir cleanup_nested threads

import argparse
import glob, os
import re
import subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import shutil


parser=argparse.ArgumentParser(description="Add full size LTR-RTs to the TE library")
parser.add_argument('EDTA_TElib_fa', type=str, help='Path to TElib')
parser.add_argument('EDTA_intact_fa', type=str, help='Path to intact lib')
parser.add_argument('work_dir', type=str, help='Path to the working directory for a sample')
parser.add_argument('cleanup_nested', type=str, help='Path to cleanup_nested.pl script')
parser.add_argument('threads', type=str, help='Number of threads to run cleanup in parallel')
args=parser.parse_args()

# Input Files
EDTA_TElib_fa = args.EDTA_TElib_fa
EDTA_intact_fa = args.EDTA_intact_fa
work_dir = args.work_dir
cleanup_nested = args.cleanup_nested
threads = args.threads

# Script dir
snake_dir = os.getcwd()

# Output Files
updated_TElib_fa = f"{os.path.splitext(EDTA_TElib_fa)[0]}.LTRupdated.fa"
final_output = f"{updated_TElib_fa}.wgn.fa"

# Genome name
gname=os.path.basename(os.path.dirname(EDTA_TElib_fa))

# Step 1: Remove LTR elements from EDTA.TElib.fa
def remove_ltr_elements(fasta_file):
    non_ltr_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if "LTR" not in record.description:
            non_ltr_records.append(record)
    return non_ltr_records

# Step 2: Cluster intact LTR elements to remove redundancy
def extract_ltr_elements(fasta_file, temp_ltr_file):
    ltr_records = []
    with open(temp_ltr_file, "w") as temp_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if "LTR" in record.description:
                ltr_records.append(record)
                SeqIO.write(record, temp_file, "fasta")
    return temp_ltr_file

def run_cleanup_nested(temp_ltr_file, threads, cleanup_nested, work_directory):
    output_file = f"{temp_ltr_file}.cln"
    os.chdir(work_directory)
    command = [
        "perl", f"{snake_dir}/{cleanup_nested}",
        "-in", temp_ltr_file,
        "-t", str(threads)
    ]
    subprocess.run(command, check=True)
    return output_file

# Step 3: Remove elements of the same family, keeping the longest
def deduplicate_families(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    family_groups = defaultdict(list)
    for record in records:
        family_id = record.id.split("|")[0]  # Adjust split logic as needed
        family_groups[family_id].append(record)
    
    representatives = []
    for family_id, group in family_groups.items():
        group.sort(key=lambda x: len(x.seq), reverse=True)
        representatives.append(group[0])  # Keep the longest
    for record in representatives:
        record.id = record.id.split("|")[0] + '#' + record.id.split("#")[1]
    filtered_representatives = [record for record in representatives if "Chr" not in record.id]
    return filtered_representatives

# Step 4: Add non-redundant LTR families to EDTA.TElib.fa
def merge_fasta(records_1, records_2):
    return records_1 + records_2

# Step 5: Update TE names
def update_te_names(records, genome_name):
    for record in records:
        record.id = re.sub(r"TE_", f"{genome_name}_TE_", record.id)
        record.description = ""  # Remove extra description if necessary
    return records

# Main Function
def main():
    # Temporary files
    temp_ltr_file = f"{os.path.splitext(EDTA_TElib_fa)[0]}_temp_LTR.intacts.fa"

    # Step 1
    non_ltr_records = remove_ltr_elements(EDTA_TElib_fa)
    
    # Step 2
    temp_ltr_file = extract_ltr_elements(EDTA_intact_fa, temp_ltr_file)
    cleaned_ltr_file = run_cleanup_nested(temp_ltr_file, threads, cleanup_nested, work_dir)
    
    # Step 3
    deduplicated_ltr_records = deduplicate_families(cleaned_ltr_file)
    
    # Step 4
    merged_records = merge_fasta(non_ltr_records, deduplicated_ltr_records)
    
    # Step 5
    genome_name = os.path.basename(EDTA_TElib_fa).split(".")[0]
    final_records = update_te_names(merged_records, genome_name)
    
    # Save Final Output
    with open(final_output, "w") as output:
        for record in final_records:
            # Write the ID and sequence without wrapping
            output.write(f">{record.id}\n{record.seq}\n")

    print(f"Updated TElib saved to {final_output}")

    # Move and rename for EDTA reannotation run
    # os.rename(EDTA_TElib_fa, f"{os.path.splitext(EDTA_TElib_fa)[0]}.originalsgEDTA.fa")
    # shutil.copy(final_output, f"{os.path.splitext(EDTA_TElib_fa)[0]}.fa")
    # shutil.copy(EDTA_TElib_fa, f"{gname}.Chr_scaffolds.fa.mod.EDTA.final")

    # Cleanup intermediate files
    os.remove(temp_ltr_file)
    os.remove(cleaned_ltr_file)
    # print(glob.glob("*TElib_temp_LTR*"))
    for f in glob.glob("*TElib_temp_LTR*"):
        os.remove(f)

# Run the Script
if __name__ == "__main__":
    main()
