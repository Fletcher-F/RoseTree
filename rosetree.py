"""Importing libraries"""
import os
import subprocess
import argparse
import Bio
from pathlib import Path
from Bio import Blast
from Bio.Seq import Seq

"""Program Info: RoseTree (Summer 2025 USRA Project)
Being developed by Fletcher Falk, this is a python program
which blasts a given input sample to NCBI and constructs a phylogeny based on results.

Metadata from NCBI is also parsed and used to supplement phylogeny output.
Started April. 28th, 2025
"""
Blast.email = "fletcherfalk@gmail.com"

"""Execute bash commands"""
def execute(command):
    result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
    print(result.stdout, "\n", result.stderr)

def main():
    """Main"""
    parser = argparse.ArgumentParser(description="RoseTree Argument List")
    parser.add_argument("--input", "-i", 
                        help="Input to be blasted against NCBI"
    )
    parser.add_argument("--threads", "-t",
                        help="Specify total number of threads to use. Default is 1.", default=1)
    args = parser.parse_args()

    """Input fasta file path to start blast"""
    if args.input:
        path = args.input
        threads= args.threads

        """Start Checking"""
        print("--- Running RoseTree v0.1 ---", "\n", "Using ", threads, " threads...", "\n", "Be sure to checkout my GitHub for other cool projects :)",
        "\n", "https://github.com/Fletcher-F", "\n",
        "File path = ", path, "\n", "\n",
        "Running qblast using Biopython on nucleotide database...."
        )

        fblast = open(path).read()
        print("Blasting input sequence against nucleotide database...", "\n")
        blast_results = Blast.qblast("blastn", "nt", fblast)

        with open("blast.xml", "wb") as out_stream:
            out_stream.write(blast_results.read())
        blast_results.close()

    """Done"""
    print("Bye!")

"""Start Program"""
if __name__ == "__main__":
    main()

"""
THINGS TO BE DONE + brainstorming:
1. Take the top results to be used in phylogeny?? Maybe based on percent match
Should I consider and how will we consider an outgroup for the tree?

2. Alignment tool with results: maybe use ClustalW or MAFFT, then put results in .phy

3. Phylogeny tool with .phy: probably IQ-tree.

4. How will we parse the metadata. Not sure if we get enough info might need to
individually search each accession number I want to use and get extra info.

5. Need to determine what metadata to add to phylogeny and how that will be done
matlab models are maybe a possibility?

6. Also need to move alot of this code into functions and clean it up.
Plus need error checking and input validation if it will be publically shared.
"""

"""Below is a bunch of extra stuff I had before for testing and stuff

Test execute, remove later
execute(f"ls {path} >> files.txt")
file_path = Path("files.txt")
file_content = file_path.read_text()
print(file_content)

file_total = len([f for f in Path(path).iterdir() if f.is_file()])

"""