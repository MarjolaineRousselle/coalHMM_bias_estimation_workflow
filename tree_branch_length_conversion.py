### This script convert a multiple newick tree file from trees with branches in number of years in trees with branches in substitution per site.
### It takes five as arguments: a mutliple newick file, a mutation rate and generation time. ( 3 arguments)
import sys
import re
from ete3 import Tree

if len(sys.argv)!=6:
	sys.exit("usage: python tree_branch_length_conversion.py <tree file> <mutation_rate> <generation_time> <tree_file_ID> <branch_ID> ")

with open(sys.argv[1]) as f:
    multiple_trees_list = f.read().splitlines() 

mutation_rate=float(sys.argv[2])
generation_time=float(sys.argv[3])
tree_file_ID=sys.argv[4]
branch_ID=sys.argv[5]

f= open('modified_tree_{}_{}.newick'.format(tree_file_ID, branch_ID),"w+")
for i in range(0, len(multiple_trees_list)):
    start=multiple_trees_list[i].partition(" ")[0]
    end=multiple_trees_list[i].partition(" ")[2].partition(" ")[0]
    real_tree=multiple_trees_list[i].partition(" ")[2].partition(" ")[2]
    t = Tree(real_tree)
    for node in t.traverse("postorder"):
      # Do some analysis on node
        node.dist=node.dist*mutation_rate*generation_time
#     print(t.write(format=1))
#     print(start)
#     print(end)
    f.write(str(start + " " + end + " " + t.write(format=1) + "\n"))


