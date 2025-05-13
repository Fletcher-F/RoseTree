"""Program Info: RoseTree (VIU Summer 2025 USRA Project)
Developed by: Fletcher Falk
April. 28th, 2025
Program blasts a given input sample to NCBI and constructs a phylogeny based on results.
Metadata from NCBI is also parsed and used to supplement phylogeny output."""

"""Importing libraries"""
import os, subprocess, argparse, sys
from pathlib import Path
from Bio import Blast, Entrez
from Bio.Seq import Seq
from roseparser import parsexml, multiparsexml, linecheck
from rosemetadata import metaparser
from rosephylogeny import phylogeny
from roseblast import nuc_blast, multi_blast

"""Function to execute bash commands"""
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
                        help="BLAST API requires an email in case of error.", required=True)
    parser.add_argument("--blasttotal", "-b",
                        help="Specify total number of blast results to use. Default is 50.", default=50)
    parser.add_argument("--multiplesequencemode", "-ms",
                        help="Turn on multiple sequence mode to allow multiple inputs. Takes a directory of fasta files and blasts all taking default results.", default=False)
    parser.add_argument("--multiblasttotal", "-mb",
                        help="Specify total number of blast results to use when running multisequencemode. Default is 10.", default=10)
    """Return arguments"""
    return parser.parse_args()

"""Runs the RoseTree program"""
def rosetree(args):
    """Assign arguments"""
    path = args.input
    threads = args.threads
    blastnum = args.blasttotal
    multiblastnum = args.multiblasttotal
    """stored as string"""
    mode = args.multiplesequencemode 
    Blast.email = args.email
    Entrez.email = args.email

    """Multiple Sequence Mode Checks"""
    if mode == "True":
        if Path(path).is_dir() == False:
            print("Error... Multiple sequence mode is enabled, therefore path requires directory.")
            sys.exit()
        files = os.listdir(path)
        for file in files:
            if os.path.splitext(file)[1] not in (".fasta", ".fa"):
                print("File in directory has invalid format...\n" "Requires .fasta or .fa")
                sys.exit()
    else: 
        """Check if file exists before starting"""
        if Path(path).is_file() == False:
            print("Error... File does not exist, exiting.")
            sys.exit()
        """Check valid fasta file path"""
        if os.path.splitext(path)[1] not in (".fasta", ".fa"):
            print("Invalid file format...\n" "Requires .fasta or .fa")
            sys.exit()

    """Start"""
    print("--- Running RoseTree v0.1 --- \n", "Using ", threads, " threads... \n", "Be sure to checkout my GitHub :) \n", 
    "https://github.com/Fletcher-F \n", "File(s) path = ", path, "\n", "Note: if rerunning program, ensure you are working in a new directory or move old output.\n",
    "Results are appended and will cause errors in downstream steps.\n\n", "Note: ensure fasta input is only 2 lines: one for title, one for sequence.\n\n"
    "Running qblast using Biopython on nucleotide database...\n")

    if mode == "True":
        print("Running in multiple sequence mode, taking top ", multiblastnum, " results from each input...\n")
        if blastnum != 50:
            print("Notice: you set the total blast number argument.\n You are running in multisequencemode",
                  " if you want to change the blast number set it instead with -mb")
        multi_blast(path, multiblastnum)
        hit_list, inputblast = multiparsexml(path)    
        """Convert input blast to string for rest of program"""
        if type(inputblast) is list:
            inputblast = ", ".join(str(x) for x in inputblast)
            inputblast = inputblast.replace(">", "")
    
    else:
        """Blast fasta to NCBI"""
        print("Taking top ", blastnum, " results.\n")
        nuc_blast(path, blastnum)
        """Parse xml into fasta"""
        hit_list, inputblast = parsexml(path)

    """Check lines written"""
    linecheck()

    """Pull metadata after parsing"""
    metaparser(hit_list, inputblast)

    """MAFFT for sequence alignment of fasta"""
    print("Performing alignment from results using MAFFT...")
    print("Using L-INS-i as less than 200 sequences are being aligned. \n")
    alignment = "mafft-linsi blastresults.fasta > alignedresults.fasta"
    execute(alignment)

    """Phylogeny construction with raxml
    Note: This tree is maximum likelihood not a bayesian tree"""
    print("Building maximum likelihood tree with RAxML all-in-one analysis...")
    mltree = "raxml-ng --all --msa alignedresults.fasta --model GTR+G --tree pars{25},rand{25} --threads " + threads
    #Long analysis option: mltree = "raxml-ng --all --msa alignedresults.fasta --model LG+G8+F --tree pars{10} --bs-trees 200 --threads " + threads
    execute(mltree)

    """Clean up files"""
    cleanup = ["rm alignedresults.fasta.raxml.r*", "rm alignedresults.fasta.raxml.mlTrees", "rm alignedresults.fasta.raxml.bootstraps", "rm alignedresults.fasta.raxml.startTree", "rm alignedresults.fasta.raxml.support" ]
    for x in cleanup:
        execute(x)

    """Draw Tree"""
    print("Drawing final tree with ETE...\n", "Note: if an error occurs you may have to manually install PyQt5: as it didn't install with ete3...\n",
    "run: pip3 install PyQt5")
    inputtree = "alignedresults.fasta.raxml.bestTree"
    phylogeny(inputtree, inputblast)

    """Done"""
    print("Final tree exported as pdf...")
    print("Bye!")

"""Main: sets up arguments and runs program"""
def main():
    args = argument_parser()
    rosetree(args)

"""Start Program"""
if __name__ == "__main__":
    main()
