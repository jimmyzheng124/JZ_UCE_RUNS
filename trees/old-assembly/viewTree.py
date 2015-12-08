import dendropy; import os, os.path; import re;
from dendropy import TaxonSet, Tree, TreeList, DnaCharacterMatrix, DataSet;
tree = dendropy.Tree.get_from_path("RAxML_cichlid.unpart.bootstrap.simoOutgroup.tre", schema="newick")
tree.print_plot();
