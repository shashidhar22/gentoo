import pandas as pd
from skbio import DistanceMatrix
from skbio.tree import nj
from ete3 import Tree, faces, TreeStyle
finch = pd.read_csv('Finch_Candida51mer_final.dist', header=0, index_col=None)
samples = finch.columns
jaccard_matrix = finch.to_numpy()
jaccard_matrix = DistanceMatrix(jaccard_matrix, samples)
newick_str = nj(jaccard_matrix, result_constructor=str)
out_file = open('Finch_Candida51mer_final.nwk', 'w')
out_file.write('{0}\n'.format(newick_str))
out_file.close()
