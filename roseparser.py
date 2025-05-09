"""XML Parser for Rosetree
Writes the results of blast from XML into FASTA for alignment and phylogeny.
Developed by: Fletcher Falk"""

from xml.etree import ElementTree as ET

"""XML Parser Function"""
def parsexml(path):
    print("Writing result sequences (100) into new fasta file for phylogeny construction...", "\n")
    xmltree = ET.parse("blast.xml")
    root = xmltree.getroot()
    hit_list = []
    
    """Read contents from input sequence to append in"""
    inputsequence = open(path, "r")
    inputblast = inputsequence.readline()

    """Open file and write input sequence"""
    newfasta = open("blastresults.fasta", "a")
    """Since readline earlier put first line into file instead of open and close"""
    newfasta.write(inputblast)
    for line in inputsequence:
        newfasta.write(line)
    newfasta.write('\n')
    inputsequence.close()

    """Write rest of blast results into file"""
    for hit in root.findall('.//Hit'):
        hit_id = hit.findtext('Hit_accession')
        hit_def = hit.findtext('Hit_def')
        hsp_seq = hit.findtext('.//Hit_hsps/Hsp/Hsp_hseq')

        """Ensure fasta format"""
        newfasta.write(">")
        newfasta.write(hit_id)
        newfasta.write(" ")
        newfasta.write(hit_def)
        newfasta.write('\n')
        newfasta.write(hsp_seq)
        newfasta.write('\n')

        """Add accession to list for metadata"""
        hit_list.append(hit_id)

    newfasta.close()

    """Check number of lines written, should be 202"""
    writecheck = open("blastresults.fasta", "r")
    checklines = len(writecheck.readlines())
    writecheck.close()

    print("Finished parsing.. ", "Wrote: ", checklines, " lines to document blastresults.fasta", "\n")
    
    """Return info for metadata after parsing"""
    return hit_list, inputblast
