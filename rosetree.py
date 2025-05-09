"""Importing libraries"""
import os, subprocess, argparse, sys, Bio
from pathlib import Path
from Bio import Blast, Entrez
from Bio.Seq import Seq
from phytreeviz import TreeViz
from roseparser import parsexml
from rosemetadata import metaparser

"""Program Info: RoseTree (VIU Summer 2025 USRA Project)
Developed by: Fletcher Falk
Program blasts a given input sample to NCBI and constructs a phylogeny based on results.
Metadata from NCBI is also parsed and used to supplement phylogeny output.

April. 28th, 2025

Function to execute bash commands"""
def execute(command):
    result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
    print(result.stdout, "\n", result.stderr)

"""Arguments for running the program"""
def argument_parser():
    """Arguments"""
    parser = argparse.ArgumentParser(description="RoseTree Argument List")
    parser.add_argument("--input", "-i", 
                        help="Input to be blasted against NCBI.", required=True)
    parser.add_argument("--threads", "-t",
                        help="Specify total number of threads to use. Default is 1.", default=1)
    parser.add_argument("--email", "-e",
                        help="BLAST API requires an email in case of error or missue.", required=True)
    return parser.parse_args()

"""Blasts fasta file to NCBI
Uses qblast from biopython library"""
def nuc_blast(path):
    blastseq = open(path).read()
    """Blast sequence to nucleotide database"""
    blast_results = Blast.qblast("blastn", "nt", blastseq, hitlist_size = 100)

    """Write results into xml"""
    with open("blast.xml", "wb") as out_stream:
        out_stream.write(blast_results.read())
    blast_results.close()

"""Runs the RoseTree program"""
def rosetree(args):
    """Assign arguments"""
    path = args.input
    threads= args.threads
    Blast.email = args.email
    Entrez.email = args.email
    
    """Check if file exists before starting"""
    validfile = Path(path)
    if validfile.is_file() == False:
        print("Error... File does not exist, exiting.")
        sys.exit()

    """Check fasta file path to start blast"""
    validfasta = os.path.splitext(args.input)
    if validfasta[1] in (".fasta", ".fa"):
        print("--- Running RoseTree v0.1 --- \n", "Using ", threads, " threads... \n", "Be sure to checkout my GitHub :) \n", 
        "https://github.com/Fletcher-F \n", "File path = ", path, "\n", "Note: if rerunning program, ensure it is in new directory or move old output as output is appended and will cause errors in downstream steps.", 
        "\n\n", "Running qblast using Biopython on nucleotide database...")

        """Blast fasta to NCBI"""
        nuc_blast(path)

        """Parse xml into fasta"""
        hit_list, inputblast = parsexml(path)

        """Pull metadata after parsing"""
        metaparser(hit_list, inputblast)

        """MAFFT for sequence alignment of fasta"""
        print("Performing alignment from results using MAFFT...")
        print("Using L-INS-i as less than 200 sequences are being aligned. \n")
        alignment = "mafft-linsi blastresults.fasta > alignedresults.fasta"
        execute(alignment)

        """Phylogeny construction
        Starting with RAxML as it is easy to run:
        Consider looking at BEAST
        Note: This tree is maximum likelihood not a bayesian tree"""
        print("Building maximum likelihood tree with RAxML all-in-one analysis...")
        mltree = "raxml-ng --all --msa alignedresults.fasta --model LG+G8+F --tree pars{10} --bs-trees 200 --threads " + threads
        execute(mltree)

        """Clean up files"""
        cleanup = ["rm alignedresults.fasta.raxml.r*", "rm alignedresults.fasta.raxml.mlTrees", "rm alignedresults.fasta.raxml.bootstraps", "rm alignedresults.fasta.raxml.startTree", "rm alignedresults.fasta.raxml.support" ]
        for x in cleanup:
            execute(x)

        """Draw Tree
        Testing with phyTreeViz tool instead of running R Studio separate
        or importing R libraries and packages"""
        print("Drawing final maximum likelihood tree with phyTreeViz...")
        inputtree = "alignedresults.fasta.raxml.bestTree"
        finaltree = TreeViz(inputtree)

        """Save final tree"""
        finaltree.savefig("FinalTree.pdf", dpi=300)
 
        """Done"""
        print("Final tree exported as pdf...")
        print("Bye!")

    else:
        """Error checking input"""
        print("Invalid file format...", "\n", "Input was", validfasta[1], "\n", "Requires .fasta or .fa")

"""Main: sets up arguments and runs program"""
def main():
    args = argument_parser()
    rosetree(args)

"""Start Program"""
if __name__ == "__main__":
    main()

"""
THINGS TO BE DONE + brainstorming:
1. Take the top results to be used in phylogeny?? Maybe based on percent match
Should I consider and how will we consider an outgroup for the tree?

2. How will we parse the metadata and what to add/how to add to phylogeny. 

3. Also need to move alot of this code into functions and clean it up.
Plus need error checking and input validation if it will be publically shared.

Extra TreeViz settings:
Excluded for now as now branch length or confidence exists for current plot
See about switch when using bayesian if we shift to that.
finaltree.show_branch_length(color="red")
finaltree.show_confidence(color="blue")
finaltree.show_scale_bar()
"""
