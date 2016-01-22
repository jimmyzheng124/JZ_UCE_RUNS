import dendropy; import os, os.path; import re; import sys; import glob;

print "Choose tree from below to view: "
print '\n'.join([os.path.basename(x) for x in glob.glob("./*malawi*")])
phyinput = raw_input("Enter tree: ")
from dendropy import TaxonSet, Tree, TreeList, DnaCharacterMatrix, DataSet;
for retry in range(0,9):
	try:
		tree = dendropy.Tree.get_from_path(phyinput, schema="newick")
		tree.print_plot()
		break
	except:
		print "Invalid tree. Try inputting your selection again:"
		phyinput = raw_input()
		pass
	if retry == 8:
		print "Let's not waste time here. Good bye."
