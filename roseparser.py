"""XML Parser for Rosetree
Writes the results of blast from XML into FASTA for alignment and phylogeny.
Developed by: Fletcher Falk"""
import os
from xml.etree import ElementTree as ET

"""Grab all blast results for Entrez"""
def blastresults(root, hit_list):
    """Write rest of blast results into file"""
    for hit in root.findall('.//Hit'):
        """Add accession to list for metadata
        Checking for dupes of accession (could happen in multi mode)"""
        if hit.findtext('Hit_accession') in hit_list:
            print ("Dupe in blast list, skipping parsing...")
            continue  
        hit_list.append(hit.findtext('Hit_accession'))
    return hit_list

"""Start new fasta file for results"""
def writefasta(path, root, hit_list):
    """Read contents from input sequence to append in"""
    inputsequence = open(path, "r")
    inputid = inputsequence.readline().replace('\n', '')
    """Open file and write input sequence"""
    newfasta = open("blastresults.fasta", "a")
    """Since readline earlier put first line into file instead of open and close"""
    newfasta.write(inputid)
    newfasta.write('\n')
    newfasta.write(inputsequence.readline())
    newfasta.write('\n')
    inputsequence.close()

    """Return final list and input sequences"""
    newfasta.close()
    hit_list = blastresults(root, hit_list)
    return hit_list, inputid

"""Write rest of blast results with associated metadata from Entrez"""
def writefinalfasta():
    xmltree = ET.parse("blastmetadata.xml")
    root = xmltree.getroot()
    finalfasta = open("blastresults.fasta", "a")
    """For each sequence in parsed xml"""
    for sequence in root.findall('.//Sequence'):
        lines = [">", sequence.findtext('Sequence_accession'), " ", sequence.findtext('Sequence_id'), "\n", sequence.findtext('full_Sequence'), "\n"]
        finalfasta.write(lines)
    finalfasta.close()

"""XML Parser Function"""
def parsexml(path):
    print("Writing blast results for Entrez...", "\n")
    xmltree = ET.parse("blast1.xml")
    root = xmltree.getroot()
    hit_list = []
    return writefasta(path, root, hit_list)

"""XML Multi Parser Function"""
def multiparsexml(path):
    print("Writing multiple sequence blast results for Entrez...", "\n")
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
