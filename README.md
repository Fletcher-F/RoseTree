# RoseTree (WIP) 🌳
## Fasta BLAST and Metadata Phylogeny Constructor Pipeline

RoseTree is a pipeline tool that will BLAST a given FASTA nucleotide input using GenBank's API,
and construct a maximum likelihood phylogenetic tree with metadata from results.

- Written in Python using Biopython and ETE3
- Parses and joins XML output into new FASTA file for alignment.
- Outputs blast results, alignment, metadata, and final phylogeny.

## Description

This is a Python project I am starting to develop as part of my NSERC Undergraduate Summer Research Award at VIU.
This is my first time working with Python and much of the project is a learning experience. There will likely be some errors and bugs that I will be working on throughout the summer.

## Upcoming Changes

- Short-hand rerun method for just the phylogeny render.
- Allow user to specify metadata to include.
- 

## Installation

Currently no package exists. To run, copy the repo and run rosetree.py locally.
Example: 
```
python3 rosetree.py -i ./test/input.fasta -t 8 -e myemail@gmail.com
```
