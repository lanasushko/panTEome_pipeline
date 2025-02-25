# conda activate tools

import argparse
from Bio import SeqIO

def main(fasta_file):
    
    # output_file=fasta_file+'.lengths'

    def is_error(command):
        try:
            command()  # Execute the command
            return False  # No error occurred
        except Exception:
            return True  # An error occurred

    sequence_lengths={}
    count=1
    for record in SeqIO.parse(fasta_file, "fasta"):
        famid=record.id.split("#")[0].strip()
        # famid=record.id.split("#")[0].strip().split('|')[0].strip()
        class_=record.id.split("#")[1].strip().split("/")[0].strip()
        if is_error(lambda: record.id.split("#")[1].strip().split("/")[1].strip()):
            superfam='NA'
        else:
            superfam=record.id.split("#")[1].strip().split("/")[1].strip()
        # if is_error(lambda: record.id.split("#")[1].strip().split("/")[2].strip()):
        #     clade=''
        # else:
        #     clade=record.id.split("#")[1].strip().split("/")[2].strip()
        sequence_lengths[record.id]=[len(record.seq),count,famid,class_,superfam]
        count+=1


    # write results
    with open(output_file, "w") as fo:
        fo.write(f"Num\tID\tLength\tfamilyid\tclass\tsupfam\n")
        for seq_id, info in sequence_lengths.items():
            fo.write(f"{info[1]}\t{seq_id}\t{info[0]}\t{info[2]}\t{info[3]}\t{info[4]}\n")

if __name__ == "__main__":
    main()