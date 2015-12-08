import dendropy; import os, os.path; import re;
from dendropy import TaxonSet, Tree, TreeList, DnaCharacterMatrix, DataSet;
tree = dendropy.Tree.get_from_path("RAxML_bestTree.best", schema="newick")
tree.print_plot();
