"""Phylogeny constructor for Rosetree
Developed by: Fletcher Falk"""
import os
from ete3 import Tree, TreeStyle, TextFace, faces, NodeStyle
from xml.etree import ElementTree as ET

"""Random color generate:
Not in this tree"""

"""Custom_layout function per ETE manual
Setup as a closure function to pass extra vars as custom_layout
seems to not like extra parameters"""
def layout_parameters(root, inputblast, samplenum):
    def custom_layout(node):
        """If node is leaf, generate TextFace of metadata that matches
        Additionally set background color to genus color generated"""
        if node.is_leaf():
            """Check if node is input sample"""
            if node.name in inputblast:
                nonlocal samplenum
                inputface = TextFace("Sample: " + (str(samplenum)), tight_text=True, fgcolor = "black", fsize = 11)
                inputface.background.color = "#BECBC1"
                faces.add_face_to_node(inputface, node, column=0, position="branch-right")
                samplenum += 1 

            """Set outgroup label"""
            if node.name == "AB021400":
                testface = TextFace("Outgroup", fsize= 7)
                faces.add_face_to_node(testface, node, column=0, position="branch-bottom")
                inputface = TextFace("Phaseolibacter flectens", tight_text=True, fgcolor = "black", fsize = 11)
                faces.add_face_to_node(inputface, node, column=0, position="branch-right")

            """Add * to collapsed nodes"""
            if node.name == "MT341881" or node.name == "MT254877" or node.name == "MT341882":
                testface = TextFace("*", fsize= 11, fgcolor = "red")
                faces.add_face_to_node(testface, node, column=0, position="branch-bottom")

            """Metadata color labels"""      
            metadata1color = {
            "Isolate-21": "#FFFC60",
            "Isolate-52": "#FFFC60",
            "Isolate-42": "#FFFC60",
            "MN577300": "white",
            "KY206826": "#FF4141",
            "MT341869": "#FFFC60",
            "JX067652": "#FFFC60",
            "JX067699": "#FFFC60",
            "NR_126305": "#FFFC60",
            "KF876191": "#FFFC60",
            "OP595649": "#FFFC60",
            "MT341882": "#FFFC60",
            "MT341881": "#FFFC60",
            "PV074719": "#4FD46A",
            "MT544516": "#4FD46A",
            "MT341868": "#4FD46A",
            "MT544538": "#4FD46A",
            "MT544528": "#4FD46A",
            "MT254877": "#4FD46A",
            "OL672320": "#4FD46A",
            "OL672322": "#4FD46A"
            }
            collapsedstars = {"OL672320", "OL672322"}
            metadata2color = {
            "Isolate-21": "#4FA3D4",
            "Isolate-52": "#4FA3D4",
            "Isolate-42": "#4FA3D4",
            "MN577300": "white",
            "KY206826": "#4FA3D4",
            "MT341869": "#4FA3D4",
            "JX067652": "#4FA3D4",
            "JX067699": "#4FA3D4",
            "NR_126305": "#4FA3D4",
            "KF876191": "#4FA3D4",
            "OP595649": "#4FA3D4",
            "MT341882": "#4FA3D4",
            "MT341881": "#4FA3D4",
            "PV074719": "#4FA3D4",
            "MT544516": "#4FA3D4",
            "MT341868": "#4FA3D4",
            "MT544538": "#4FA3D4",
            "MT544528": "#4FA3D4",
            "MT254877": "#4FA3D4",
            "OL672320": "#4FA3D4",
            "OL672322": "#4FA3D4"
            }
            if node.name in metadata1color:
                if node.name in collapsedstars:
                    faces.add_face_to_node(TextFace(" *", fsize=11), node, column=1, aligned=True)
                    faces.add_face_to_node(TextFace(" *", fsize=11), node, column=3, aligned=True)
                else:
                    spacer = TextFace(" ", fsize=11)
                    spacer.background.color = "white"
                    spacer.margin_left = 5
                    spacer.margin_right = 5
                    faces.add_face_to_node(spacer, node, column=1, aligned=True)
                    faces.add_face_to_node(spacer, node, column=3, aligned=True)

                face2 = TextFace("      ", fsize=11)
                face2.background.color = metadata1color[node.name]
                face2.margin_left = 13
                face2.margin_right = 13
                face2.border.width = 1
                faces.add_face_to_node(face2, node, column=2, aligned=True)

                face3 = TextFace("      ", fsize=11)
                face3.background.color = metadata2color[node.name]
                face3.margin_left = 13
                face3.margin_right = 13
                face3.border.width = 1
                faces.add_face_to_node(face3, node, column=4, aligned=True)
                
            """For each sequence in parsed xml"""
            """Old metadata stuff"""
            for sequence in root.findall('.//Sequence'):
                if node.name == sequence.findtext('Sequence_accession'):
                    face = TextFace(sequence.findtext('Sequence_id'), tight_text=True, fsize = 11)
                    faces.add_face_to_node(face, node, column=0, position="branch-right")

    return custom_layout

"""Phylogeny Function"""
def main():
    print("Drawing final maximum likelihood tree...\n")
    """If adding metadata include this .xml
    See rosephylogeny for metadata example"""
    xmltree = ET.parse("blastmetadata.xml")
    root = xmltree.getroot()

    """Open RaXML tree file and style tree"""
    with open ("alignedresults.fasta.raxml.supportTBE", "r") as file:
        s = file.read()
        t = Tree(s)
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.show_branch_length = False
        ts.branch_vertical_margin = 1
        ts.margin_left = 25
        ts.margin_right = 25
        ts.margin_top = 25
        ts.margin_bottom = 50

        """For each leaf add leaf name manually so it can be offset color"""
        for leaf in t.iter_leaves():
            face = TextFace(leaf.name, fsize=10, fgcolor="#333333")
            leaf.add_face(face, column=0, position="branch-right")

        """For each node in the tree set node style and add support and distance values"""
        for node in t.traverse():
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "black"
            nstyle["vt_line_width"] = 1
            nstyle["hz_line_width"] = 1
            nstyle["hz_line_color"] = "black"
            nstyle["vt_line_color"] = "black"
            node.set_style(nstyle)

            """If distance < 0.0001 dont include distance value"""
            if node.dist < 0.0001:
                node.add_face(TextFace(("      "), fsize = 7, fgcolor="black"), column=0, position="branch-top")
            else:
                node.add_face(TextFace(str(node.dist), fsize = 7, fgcolor="black"), column=0, position="branch-top")

            """Adding support values
            TBE is in percent so convert to whole number"""
            if not node.is_leaf():
                if (node.support)*100 >= 60: 
                    node.add_face(TextFace(str(round((node.support*100), 0)), fsize=7, fgcolor="red"), column=0, position="branch-bottom")

            """Nodes to collapse from tree"""
            collapsenodes = ("MT341880", "OP595686", "MT341876", "MT254887", "MT254874", "MT254860", "MT341867", "MT341878")
            if node.name in collapsenodes:
                node.delete()
    
        os.environ["QT_QPA_PLATFORM"] = "offscreen"
        ts.layout_fn = layout_parameters(root, ("Isolate-21", "Isolate-52", "Isolate-42"), 1)
        t.render("outgroup.pdf", tree_style = ts)
        file.close()

"""Start Program"""
if __name__ == "__main__":
    main()

