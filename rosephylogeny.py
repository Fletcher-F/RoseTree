"""Phylogeny constructor for Rosetree
Uses ETE3 to build a phylogeny and annotate with parsed metadata
Developed by: Fletcher Falk"""
import os
from random import randint
from ete3 import Tree, TreeStyle, TextFace, faces
from xml.etree import ElementTree as ET

"""Random color generate"""
def graphcolor(colortotal):
    colorlist = []
    for x in colortotal:
        """Generates random light tone color for legibility"""
        colorlist.append(('#%02x%02x%02x' % (randint(197, 240), randint(197, 240), randint(197, 240))).upper())
    return colorlist

"""Species list and color"""
def specieslist(xmltree):
    root = xmltree.getroot()
    specieslist = []
    """For each species in blast generate list of unique genus"""
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

"""Custom_layout function per ETE manual
Setup as a closure function to pass extra vars as custom_layout
seems to not like extra parameters"""
def layout_parameters(root, colors, cleanlist):
    def custom_layout(node):
        """If node is leaf, generate TextFace of metadata that matches
        Additionally set background color to genus color generated"""
        if node.is_leaf():
            for sequence in root.findall('.//Sequence'):
                if node.name == sequence.findtext('Sequence_accession'):
                    id = sequence.findtext('Sequence_id').split(" ")
                    if id[0] == "uncultured":
                        position = cleanlist.index(id[1])
                    else:
                        position = cleanlist.index(id[0])
                    face = TextFace(sequence.findtext('Sequence_id'), tight_text=True, fsize = 12)
                    face.background.color = colors[position]
                    faces.add_face_to_node(face, node, column=0, position="branch-right")
                    face2 = TextFace(sequence.findtext('.//Source/isolation_source'), fsize = 10)
                    faces.add_face_to_node(face2, node, column=1, aligned = True)
                    face3 = TextFace(sequence.findtext('.//Source/geo_loc_name'), fsize = 10)
                    faces.add_face_to_node(face3, node, column=2, aligned = True)
                    """
                    Additional metadata that is removed for now.
                    face4 = TextFace(sequence.findtext('.//Source/db_xref'), fsize = 10)
                    faces.add_face_to_node(face4, node, column=3, aligned = True)
                    face5 = TextFace(sequence.findtext('.//Source/mol_type'), fsize = 10)
                    faces.add_face_to_node(face5, node, column=4, aligned = True)
                    """
    return custom_layout

"""Phylogeny Function"""
def phylogeny(inputtree, inputblast):
    print("Drawing final maximum likelihood tree with ETE...\n")

    """Construct variables for custom_layout"""
    xmltree = ET.parse("blastmetadata.xml")
    root = xmltree.getroot()
    colors, cleanlist = specieslist(xmltree)

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
        styletree.title.add_face(TextFace("Constructed Maximum Likelihood Phylogeny from NCBI Blast of: " + inputblast, fsize=16), column=0)
        """Render tree
        Note: os.environ is used due to render bug"""
        os.environ["QT_QPA_PLATFORM"] = "offscreen"
        styletree.layout_fn = layout_parameters(root, colors, cleanlist)
        tree.render("final-tree-render.pdf", tree_style = styletree)
    treefile.close()
