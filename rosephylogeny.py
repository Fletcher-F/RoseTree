"""Phylogeny constructor for Rosetree
Uses ETE3 to build a phylogeny and annotate with parsed metadata
Developed by: Fletcher Falk"""
import os, re
from random import randint
from ete3 import Tree, TreeStyle, TextFace, faces
from xml.etree import ElementTree as ET

"""Random color generate"""
def graphcolor(colortotal):
    colorlist = []
    for x in colortotal:
        """Generates random calm/lighter tone color for legibility"""
        colorlist.append(('#%02x%02x%02x' % (randint(197, 240), randint(197, 240), randint(197, 240))).upper())
    return colorlist

"""Species list and color"""
def specieslist(xmltree):
    root = xmltree.getroot()
    specieslist = []
    """For each species in blast generate list of unique genus to generate colors for"""
    for sequence in root.findall('.//Sequence'):
        id = sequence.findtext('Sequence_id')
        tmplist = id.split(" ")
        if tmplist[0] == "uncultured":
            specieslist.append(tmplist[1])
        else:
            specieslist.append(tmplist[0])
    cleanlist = list(set(specieslist))
    colors = graphcolor(cleanlist)
    """Return list and color for each genus"""
    return colors, cleanlist

"""Node text face fixer"""
def textfacefix(inputlabel):
    """If input is longer than 25 characters add newline to format better"""
    max_length = 25
    finalabel = ""
    """When inputlabel is greater than 25 split and reformat"""
    while len(inputlabel) > max_length:
        whitespace = re.search(r"\s", inputlabel[max_length:])
        """If whitespace is found find split point otherwise use string"""
        if whitespace:
            split = max_length + whitespace.start()
        else:
            split = len(inputlabel)
        """Remove whitespace and add newline"""
        finalabel += inputlabel[:split].strip() + "\n"
        inputlabel = inputlabel[split:].strip()
    finalabel += inputlabel
    return finalabel

"""Custom_layout function per ETE manual
Setup as a closure function to pass extra vars as custom_layout
seems to not like extra parameters"""
def layout_parameters(root, colors, cleanlist, inputblast, samplenum):
    def custom_layout(node):
        """If node is leaf, generate TextFace of metadata that matches
        Additionally set background color to genus color generated"""
        if node.is_leaf():
            """Check if node is one of the input fasta files"""
            if node.name in inputblast:
                nonlocal samplenum
                inputface = TextFace("Sample: " + (str(samplenum)), tight_text=True, fgcolor = "white", fsize = 12)
                inputface.background.color = "#F7879A"
                faces.add_face_to_node(inputface, node, column=0, position="branch-right")
                samplenum += 1 
            """For each sequence in parsed xml"""
            for sequence in root.findall('.//Sequence'):
                if node.name == sequence.findtext('Sequence_accession'):
                    id = sequence.findtext('Sequence_id').split(" ")
                    """Get color for node"""
                    if id[0] == "uncultured":
                        position = cleanlist.index(id[1])
                    else:
                        position = cleanlist.index(id[0])
            
                    """Add text faces for each metadata label"""
                    face = TextFace(sequence.findtext('Sequence_id'), tight_text=True, fsize = 12)
                    face.background.color = colors[position]
                    faces.add_face_to_node(face, node, column=0, position="branch-right")
                    face2 = TextFace(textfacefix(str(sequence.findtext('.//Source/isolation_source'))), fgcolor = "gray", fsize = 10)
                    face2.margin_left = 15
                    face2.margin_bottom = 10
                    faces.add_face_to_node(face2, node, column=1, aligned = True)
                    face3 = TextFace(textfacefix(str(sequence.findtext('.//Source/geo_loc_name'))), fgcolor = "gray", fsize = 10)
                    face3.margin_left = 15
                    face3.margin_bottom = 10
                    faces.add_face_to_node(face3, node, column=2, aligned = True)
                    face4 = TextFace(textfacefix(str(sequence.findtext('.//Source/host'))), fgcolor = "gray", fsize = 10)
                    face4.margin_left = 15
                    face4.margin_bottom = 10
                    faces.add_face_to_node(face4, node, column=3, aligned = True)
                    face5 = TextFace(textfacefix(str(sequence.findtext('.//Source/db_xref'))), fgcolor = "gray", fsize = 10)
                    face5.margin_left = 15
                    face5.margin_bottom = 10
                    faces.add_face_to_node(face5, node, column=4, aligned = True)
    return custom_layout

"""Phylogeny Function"""
def phylogeny(inputtree, inputblast):
    print("Drawing final maximum likelihood tree with ETE...\n")

    """Construct variables for custom_layout"""
    xmltree = ET.parse("blastmetadata.xml")
    root = xmltree.getroot()
    colors, cleanlist = specieslist(xmltree)
    samplenum = 1

    """Open best tree from raxml"""
    with open (inputtree, "r") as treefile:
        tree = Tree(treefile.read())
        """Tree Styling"""
        styletree = TreeStyle()
        styletree.show_leaf_name = True
        styletree.show_branch_length = True
        styletree.extra_branch_line_type = 0
        styletree.margin_left = 30
        styletree.margin_right = 30
        styletree.margin_top = 30
        styletree.margin_bottom = 30

        """Tree Title"""
        styletree.title.add_face(TextFace("Maximum Likelihood Phylogeny from NCBI Blast of:\n" + inputblast, fsize=14), column=0)
        """Render tree
        Note: os.environ is used due to render bug"""
        os.environ["QT_QPA_PLATFORM"] = "offscreen"
        styletree.layout_fn = layout_parameters(root, colors, cleanlist, inputblast, samplenum)
        tree.render("final-tree-render.pdf", tree_style = styletree)
    treefile.close()
