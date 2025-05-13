"""XML Parser for Rosetree
Writes the results of blast from XML into FASTA for alignment and phylogeny.
Developed by: Fletcher Falk"""

import os
from xml.etree import ElementTree as ET

"""Write new fasta file parsed from XML"""
def writefasta(path, root, hit_list):
    """Read contents from input sequence to append in"""
    inputsequence = open(path, "r")
    inputblast = inputsequence.readline()

    """Open file and write input sequence"""
    newfasta = open("blastresults.fasta", "a")
    """Since readline earlier put first line into file instead of open and close"""
    newfasta.write(inputblast)
    newfasta.write(inputsequence.readline())
    newfasta.write('\n')
    inputsequence.close()

    """Write rest of blast results into file"""
    for hit in root.findall('.//Hit'):

        """Add accession to list for metadata
        Checking for dupes of accession (could happen in multi mode)"""
        if hit.findtext('Hit_accession') in hit_list:
            print ("Dupe in blast list, skipping parsing...")
            continue
        
        hit_list.append(hit.findtext('Hit_accession'))
        """Ensure fasta format"""
        newfasta.write(">")
        newfasta.write(hit.findtext('Hit_accession'))
        newfasta.write(" ")
        newfasta.write(hit.findtext('Hit_def'))
        newfasta.write('\n')
        newfasta.write(hit.findtext('.//Hit_hsps/Hsp/Hsp_hseq'))
        newfasta.write('\n')
    """Return final list and input sequences"""
    newfasta.close()
    return hit_list, inputblast

"""XML Parser Function"""
def parsexml(path):
    print("Writing result sequences into new fasta file for phylogeny construction...", "\n")
    xmltree = ET.parse("blast.xml")
    root = xmltree.getroot()
    hit_list = []
    return writefasta(path, root, hit_list)

"""XML Multi Parser Function"""
def multiparsexml(path):
    print("Writing multiple sequence results into new fasta file for phylogeny construction...", "\n")
    blastseqs = [os.path.join(path, file) for file in os.listdir(path)]
    hit_list = []
    inputsequences = []
    blastnum = 1
    """For each blast sequence parse info into xml"""
    for blast in blastseqs:
        xmltree = ET.parse(("blast" + str(blastnum) + ".xml"))
        root = xmltree.getroot()
        sequence = writefasta(blast, root, hit_list)[1]
        inputsequences.append(sequence)
        blastnum += 1
    """Return final list and input sequences"""
    return hit_list, inputsequences

def linecheck():
    """Check number of lines written, should be (2*blast number) + 2"""
    writecheck = open("blastresults.fasta", "r")
    checklines = len(writecheck.readlines())
    writecheck.close()

    print("Finished parsing.. ", "Wrote: ", checklines, " lines to document blastresults.fasta", "\n")  
    """Return info for metadata after parsing"""
