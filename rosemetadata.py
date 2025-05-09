"""Metadata Parser for Rosetree
Uses Entrez to pull associated data for each accession number from BLAST.
Developed by: Fletcher Falk"""

import Bio, time
from Bio import Entrez
from xml.etree import ElementTree as ET

"""Entrez Metadata Function"""
def metaparser(accessionlist, inputblast):
    print("Fetching metadata with Entrez from BLAST results...", "\n")

    """Start writing new xml output for metadata by building ElementTree"""
    root = ET.Element("MetadataOutput")
    ET.SubElement(root, "Info",).text = "Rosetree metadata from blast of " + inputblast
    ET.SubElement(root, "Dev",).text = "Part of program made by Fletcher Falk"
    """For sequence number"""
    count = 1

    """Go through each accession from the input list.
    Setup like this to allow for sleep and end of loop to follow NCBI Guidelines
    and prevent overloading requests of all the list at once."""
    for accession in accessionlist:
        """Fetch with Entrez"""
        fetch = Entrez.efetch(db = "nucleotide", id = accession, rettype = "gb", retmode = "xml")
        metadata = Entrez.read(fetch)
        fetch.close()

        sequence = ET.SubElement(root, "Sequence")

        """Add a new element off number"""
        ET.SubElement(sequence, "Sequence_num").text = str(count)
        ET.SubElement(sequence, "Sequence_id").text = metadata[0]['GBSeq_organism']
        """Add source details based on what is available for the given genbank accession"""
        for feature in metadata[0]['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'source':
                source = ET.SubElement(sequence, "Source")
                for qual in feature['GBFeature_quals']:
                    if qual['GBQualifier_name'] == 'mol_type':
                        ET.SubElement(source, "mol_type").text = qual['GBQualifier_value']
                    if qual['GBQualifier_name'] == 'isolation_source':
                        ET.SubElement(source, "isolation_source").text = qual['GBQualifier_value']
                    if qual['GBQualifier_name'] == 'geo_loc_name':
                        ET.SubElement(source, "geo_loc_name").text = qual['GBQualifier_value']
                    if qual['GBQualifier_name'] == 'db_xref':
                        ET.SubElement(source, "db_xref").text = qual['GBQualifier_value']
        count += 1

    """Sleep to follow appropriate NCBI Guidelines for requests"""
    time.sleep(0.5)

    """Output final xml file"""
    outputxml = ET.ElementTree(root)
    outputxml.write("blastmetadata.xml")