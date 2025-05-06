"""Importing libraries"""
import os
import subprocess
import argparse
import Bio
from pathlib import Path
from Bio import Blast
from Bio.Seq import Seq
from xml.etree import ElementTree as ET

"""Program Info: RoseTree (Summer 2025 USRA Project)
Being developed by Fletcher Falk, this is a python program
which blasts a given input sample to NCBI and constructs a phylogeny based on results.

Metadata from NCBI is also parsed and used to supplement phylogeny output.
Started April. 28th, 2025
"""

"""Change blast email to a argument later"""
Blast.email = "fletcherfalk@gmail.com"

"""Execute bash commands"""
def execute(command):
    result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
    print(result.stdout, "\n", result.stderr)

"""XML Parser function for pulling results from blast"""
def parsexml(path):
    print("Writing top 30 results into new fasta file for phylogeny construction...", "\n")
    xmltree = ET.parse("blast.xml")
    root = xmltree.getroot()

    accessionList = []
    titleList = []
    seqList = []

    """I would like to do this without 3 for loops but it was being weird"""
    for Hit_accession in root.iter('Hit_accession'):
        accessionList.append(Hit_accession.text)

    for Hit_def in root.iter('Hit_def'):
        titleList.append(Hit_def.text)

    for Hsp_hseq in root.iter('Hsp_hseq'):
        seqList.append(Hsp_hseq.text)

    """Read contents from input sequence to append in"""
    inputsequence = open(path, "r")

    """Open file and write in input sequence"""
    newfasta = open("blastresults.fasta", "a")
    for line in inputsequence:
        newfasta.write(line)
    newfasta.write('\n')

    """Write top 30 blast results"""
    for x in range(30):
        newfasta.write(">")
        newfasta.write(accessionList[x])
        newfasta.write(titleList[x])
        newfasta.write('\n')
        newfasta.write(seqList[x])
        newfasta.write('\n')
    newfasta.close()
    """Close, done"""
    writecheck = open("blastresults.fasta", "r")
    checklines = len(writecheck.readlines())
    writecheck.close()

    print("Finished parsing.. ", "Wrote: ", checklines, " lines to document blastresults.fasta", "\n", "\n")

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
        "Running qblast using Biopython on nucleotide database..."
        )

        fblast = open(path).read()
        """Blast sequence to nucleotide database"""
        blast_results = Blast.qblast("blastn", "nt", fblast)

        """Write results into xml"""
        with open("blast.xml", "wb") as out_stream:
            out_stream.write(blast_results.read())
        blast_results.close()

        """Parse xml into fasta"""
        parsexml(path)

        """MAFFT for sequence alignment of fasta"""
        print("Performing alignment from results using MAFFT...", "\n")
        print("Using L-INS-i as less than 200 sequences are being aligned.", "\n")
        alignment = "mafft-linsi blastresults.fasta > alignedresults.fasta"
        execute(alignment)

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
