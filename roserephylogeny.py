"""Phylogeny constructor for Rosetree
Uses ETE3 to build a phylogeny and annotate with parsed metadata
This can be ran after the complete pipeline if you want to make edits and rerun the
rendering.
Developed by: Fletcher Falk"""
import os, argparse
from rosephylogeny import phylogeny

"""Arguments for running the program"""
def argument_parser():
    """Arguments"""
    parser = argparse.ArgumentParser(description="RosePhylogeny Argument List")
    parser.add_argument("--input", "-i", 
                        help="Rerender the phylogeny using the input tree file", required=True)
    """Return arguments"""
    return parser.parse_args()

def roserephylogeny(args):
    """Assign arguments"""
    input = (args.input)
    """Draw Tree"""
    print("Drawing final tree with ETE...\n", "Note: if an error occurs you may have to manually install PyQt5: as it didn't install with ete3...\n",
    "run: pip3 install PyQt5")

    """Phylogeny data for rerun"""
    phydata = open("phydata", "r")
    inputblast = phydata.readline()
    phylogeny(input, inputblast, "final-tree-rerender")

    """Done"""
    print("Final tree exported as pdf...")
    print("Bye!")

"""Main: sets up arguments and runs program"""
def main():
    args = argument_parser()
    roserephylogeny(args)

"""Start Program"""
if __name__ == "__main__":
    main()