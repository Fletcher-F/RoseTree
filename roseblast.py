"""NCBI Blast for Rosetree
Writes the results of blast into XML for parsing.
Developed by: Fletcher Falk"""

import os, time
from Bio import Blast

"""Blasts fasta file to NCBI
Uses qblast from Biopython library"""
def nuc_blast(path, blastnum):
    blastseq = open(path).read()
    """Blast sequence to nucleotide database"""
    blast_results = Blast.qblast("blastn", "nt", blastseq, hitlist_size = blastnum)

    """Write results into xml"""
    with open("blast.xml", "wb") as out_stream:
        out_stream.write(blast_results.read())
    blast_results.close()

"""Blasts fasta files to NCBI
takes directory path instead"""
def multi_blast(path, blastnum):
    blastseqs = [os.path.join(path, file) for file in os.listdir(path)]
    blastloop = 1
    for fasta in blastseqs:
        """Blast each fasta to NCBI"""
        blastseq = open(fasta).read()
        blast_results = Blast.qblast("blastn", "nt", blastseq, hitlist_size = blastnum)
        
        """Write results into xml"""
        with open("blast" + str(blastloop) + ".xml", "wb") as out_stream:
            out_stream.write(blast_results.read())
        blast_results.close()
        blastloop += 1
        """Fit NCBI API guidelines"""
        time.sleep(10)