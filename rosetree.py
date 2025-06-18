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
from roseparser import parsexml, multiparsexml, linecheck, writefinalfasta
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
    parser.add_argument("--outgroup", "-o",
                        help="Specify the outgroup to be used for the phylogeny. (Give it a path to a fasta file with sequence for outgroup)")
    parser.add_argument("--blasttotal", "-b",
                        help="Specify total number of blast results to use. Default is 50.", default=50)
    parser.add_argument("--multiplesequencemode", "-ms",
                        help="Turn on multiple sequence mode (True) to allow multiple inputs. Takes a directory of fasta files and blasts all taking default results.", default=False)
    parser.add_argument("--multiblasttotal", "-mb",
                        help="Specify total number of blast results to use when running multisequencemode. Default is 10.", default=10)
    """Return arguments"""
    return parser.parse_args()

"""Check valid fasta"""
def valid_fasta(file):
    if os.path.splitext(file)[1] not in (".fasta", ".fa"):
        print("Invalid file format...\n" "Requires .fasta or .fa")
        sys.exit()

"""Single blast mode"""
def single_mode(path, blastnum, multiblastnum = 10):
    """Check if file exists before starting"""
    if Path(path).is_file() == False:
        print("Error... File does not exist, exiting.")
        sys.exit()
    """Check valid fasta file path"""
    valid_fasta(path)

    """Blast single fasta to NCBI"""
    print("Taking top ", blastnum, " results.\n")
    nuc_blast(path, blastnum, 1)
    """Parse xml into fasta"""
    return parsexml(path)

"""Multi blast mode"""
def multi_mode(path, blastnum, multiblastnum):
    """Check if directory exists before starting"""
    if Path(path).is_dir() == False:
        print("Error... Multiple sequence mode is enabled, therefore path requires directory.")
        sys.exit()
    files = os.listdir(path)
    """Check valid fasta files in directory"""
    for file in files:
        valid_fasta(file)
    
    """Blast multi fasta to NCBI"""
    print("Running in multiple sequence mode, taking top ", multiblastnum, " results from each input...\n")
    if blastnum != 50:
        print("Notice: you set the total blast number argument.\n You are running in multisequencemode",
        " if you want to change the blast number set it instead with -mb")
    multi_blast(path, multiblastnum)
    return multiparsexml(path)

"""Runs the RoseTree program"""
def rosetree(args):
    """Assign arguments"""
    path = args.input
    threads = args.threads
    blastnum = args.blasttotal
    multiblastnum = args.multiblasttotal
    outgroup = args.outgroup
    """stored as string"""
    mode = args.multiplesequencemode 
    Blast.email = args.email
    Entrez.email = args.email

    """Start"""
    print("--- Running RoseTree v0.1 --- \n", "Using ", threads, " threads... \n", "Be sure to checkout my GitHub :) \n", 
    "https://github.com/Fletcher-F \n", "File(s) path = ", path, "\n", "Note: if rerunning program, ensure you are working in a new directory or move old output.\n",
    "Results are appended and will cause errors in downstream steps.\n\n", "Note: ensure fasta input is only 2 lines: one for title, one for sequence.\n\n"
    "Running qblast using Biopython on nucleotide database...\n")

    """Checking outgroup input
    Then write to fasta for alignment"""
    if outgroup:
        print("You specified an outgroup, checking valid format..\n")
        valid_fasta(outgroup)
        outgroup = open(outgroup, "r")
        initialfasta = open("blastresults.fasta", "a")
        outgroupname = outgroup.readline()
        initialfasta.write(outgroupname)
        outgroupname = "--outgroup " + outgroupname.split(" ")[0].replace(">", "")
        initialfasta.write(outgroup.readline())
        initialfasta.write("\n")
        initialfasta.close()
    else:
        print("WARNING: No outgroup was given. Recommend restarting the program with an outgroup (-o ./test/outgroup.fasta\n)")
        outgroupname = None

    modes = {
        "False": single_mode,
        "True": multi_mode
    }
    hit_list, inputblast = modes.get(mode, single_mode)(path, blastnum, multiblastnum)

    """Convert input blast to string for downstream use
    returned as list if given multiple fastas in multiple sequence mode"""
    if type(inputblast) is list:
        inputblast = (", ".join(str(x) for x in inputblast))
    inputblast = inputblast.replace(">", "")
    
    """Pull metadata after parsing"""
    metaparser(hit_list, inputblast)

    """Write rest of blast results into fasta after pulling complete sequence with Entrez"""
    writefinalfasta()

    """Check lines written"""
    linecheck()

    """MAFFT for sequence alignment of fasta"""
    print("Performing alignment from results using MAFFT...")
    alignment = "mafft --auto --thread %s blastresults.fasta > alignedresults.fasta" % threads
    execute(alignment)

    """Trim alignment with Trimal
    Currently disabled: trimming is a bit too harsh for my current analysis
    print("Trimming alignment results using Trimal...")
    Automated1 option is optimized for ML trees
    trim = "trimal -in alignedresults.fasta -out trimmedalignedresults.fasta -automated1"
    execute(trim)
    """

    """Perform model test on alignment"""
    print("Running model test on alignment results...")
    modeltest = "modeltest-ng -i alignedresults.fasta -t ml -p %s -r 12345" % threads
    execute(modeltest)

    optimalmodel = open("alignedresults.fasta.out", "r")
    for line in optimalmodel:
        if "raxml-ng" in line:
            commands = line.split(" ")
            command = commands[-1].rstrip("\n")
            print("Model selected: ", command)
            break

    """Phylogeny construction with RaXML"""
    print("Building maximum likelihood tree with RAxML all-in-one analysis...")
    mltree = "raxml-ng --all --msa alignedresults.fasta --model %s --tree pars{25},rand{25} --bs-trees 2500 --bs-metric fbp,tbe --seed 12345 --threads %s %s" % (command, threads, outgroupname)
    execute(mltree)

    """Draw Tree"""
    print("Drawing final tree (FBP support) with ETE...\n", "Note: if an error occurs you may have to manually install PyQt5: as it didn't install with ete3...\n",
    "run: pip3 install PyQt5")
    inputtree = "alignedresults.fasta.raxml.supportFBP"
    """Phylogeny data for rerun"""
    phydata = open("phydata", "a")
    phydata.write(inputblast)
    phylogeny(inputtree, inputblast, "final-tree-render")

    """Done"""
    print("Final tree exported as pdf...")
    print("Check the model test results and confirm model used.")
    print("If you would like to rerender, use rephylogeny.")
    print("Bye!")

"""Main: sets up arguments and runs program"""
def main():
    args = argument_parser()
    rosetree(args)

"""Start Program"""
if __name__ == "__main__":
    main()
