"""Metadata Parser for Rosetree
Uses Entrez to pull associated data for each accession number from BLAST.

!!! Alternative version that removes duplicates that have same metadata, see below for more info. !!!

Developed by: Fletcher Falk"""
import Bio, time
from Bio import Entrez
from xml.etree import ElementTree as ET

"""Entrez formatting 16S from complete genome sequence"""
def parse16S(metadata, accession):
    good16s = False
    sequence = ""
    """Look through feature table until 16S rRNA is found"""
    for feature in metadata[0]['GBSeq_feature-table']:
        if sequence != "":
            break
        if feature['GBFeature_key'] == 'rRNA':
            for qual in feature['GBFeature_quals']:
                if qual['GBQualifier_name'] == "product" and "16S ribosomal RNA" in qual['GBQualifier_value']:
                    good16s = True
                elif qual['GBQualifier_name'] == "transcription" and good16s == True:
                    sequence = qual['GBQualifier_value']
                    break
    if sequence == "":
        print("Unable to grab sequence from complete genome ", accession, " Likely unannotated.. excluding") 
        return "Error"
    return sequence

"""Entrez Metadata Function"""
def metaparser(accessionlist, inputblast):
    print("Fetching metadata with Entrez from BLAST results...", "\n")

    """Start writing new xml output for metadata by building ElementTree"""
    root = ET.Element("MetadataOutput")
    ET.SubElement(root, "Info",).text = "Rosetree metadata from blast of " + inputblast
    ET.SubElement(root, "Dev",).text = "Part of program made by Fletcher Falk"
    """For sequence number"""
    count = 1

    duplist = []

    """Go through each accession from the input list.
    Setup like this to allow for sleep at end of loop to follow NCBI Guidelines
    and prevent overloading requests of all the list at once."""
    for accession in accessionlist:
        """Fetch with Entrez"""
        fetch = Entrez.efetch(db = "nucleotide", id = accession, rettype = "gb", retmode = "xml")
        metadata = Entrez.read(fetch)
        fetch.close()

        host = "blank"
        isosource = "blank"
        authors = "Unknown"
        
        """REMOVING DUPLICATE SEQUENCES"""
        for feature in metadata[0]['GBSeq_references']:
            if feature['GBReference_reference'] == '1':
                authors = str(feature.get('GBReference_authors'))

        for feature in metadata[0]['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'source':
                for qual in feature['GBFeature_quals']:
                    match qual['GBQualifier_name']:
                        case 'isolation_source':
                            isosource = qual['GBQualifier_value']
                        case 'host':
                            host = qual['GBQualifier_value']
        metalist = authors+" "+isosource+" "+host

        if str(metalist) in duplist:
            print("Skipping dupe in list...")
            continue
        else:
            duplist.append(str(metalist))

            """Add full sequence for phylogeny"""
            if "complete genome" or "chromosome" in str(metadata[0]['GBSeq_definition']):
                rRNAsequence = parse16S(metadata, accession)
                if rRNAsequence == "Error":
                    continue
                sequence = ET.SubElement(root, "Sequence")
                ET.SubElement(sequence, "full_Sequence").text = rRNAsequence
            else:
                sequence = ET.SubElement(root, "Sequence")
                ET.SubElement(sequence, "full_Sequence").text = metadata[0]['GBSeq_sequence']
            
            """Add a new element off number"""
            ET.SubElement(sequence, "Sequence_num").text = str(count)
            ET.SubElement(sequence, "Sequence_accession").text = accession
            ET.SubElement(sequence, "Sequence_id").text = metadata[0]['GBSeq_organism']

            for feature in metadata[0]['GBSeq_references']:
                if feature['GBReference_reference'] == '1':
                    ET.SubElement(sequence, "Authors").text = str(feature.get('GBReference_authors'))

            """Add source details based on what is available for the given genbank accession"""
            for feature in metadata[0]['GBSeq_feature-table']:
                if feature['GBFeature_key'] == 'source':
                    source = ET.SubElement(sequence, "Source")
                    for qual in feature['GBFeature_quals']:
                        match qual['GBQualifier_name']:
                            case 'mol_type':
                                ET.SubElement(source, "mol_type").text = qual['GBQualifier_value']
                            case 'isolation_source':
                                ET.SubElement(source, "isolation_source").text = qual['GBQualifier_value']
                            case 'host':
                                ET.SubElement(source, "host").text = qual['GBQualifier_value']
                            case 'geo_loc_name':
                                ET.SubElement(source, "geo_loc_name").text = qual['GBQualifier_value']
                            case 'db_xref':
                                ET.SubElement(source, "db_xref").text = qual['GBQualifier_value']
            count += 1
        """Sleep to follow appropriate NCBI Guidelines for requests"""
        time.sleep(0.5)

    """Output final xml file"""
    outputxml = ET.ElementTree(root)
    outputxml.write("blastmetadata.xml")
