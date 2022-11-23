###########################################################################################
## This script contains the code to reproduce plots in Ext Fig 3B                        ##
## Mensah & Niskanen et al.                                                              ##
## Aberrant phase separation and nucleolar dysfunction in human genetic disease 2022     ##
## Author: Alexandre P Magalhaes                                                         ##
###########################################################################################

import warnings
import re
import sys
import glob
import numpy as np
import pickle
from alphafold.data import pipeline

import py3Dmol
import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42

args = sys.argv
dir = args[1]
pkls = glob.glob(dir + '/*model*.pkl')
msa_file = dir + '/msas/bfd_uniclust_hits.a3m'

homooligomer = 0

row = 0
with open(msa_file, 'r') as m:
	for line in m:
		row += 1
		if  row == 2:
			query_sequence = line.strip()
			break


a3m_lines = "".join(open(msa_file,"r").readlines())
msa, deletion_matrix = pipeline.parsers.parse_a3m(a3m_lines)

deduped_full_msa = list(dict.fromkeys(msa))
msa_arr = np.array([list(seq) for seq in deduped_full_msa])
seqid = (np.array(list(query_sequence)) == msa_arr).mean(-1)
seqid_sort = seqid.argsort() #[::-1]
non_gaps = (msa_arr != "-").astype(float)
non_gaps[non_gaps == 0] = np.nan

plt.figure(figsize=(14,4),dpi=300)

plt.subplot(1,2,1); plt.title("Sequence coverage")
plt.imshow(non_gaps[seqid_sort]*seqid[seqid_sort,None],
           interpolation='nearest', aspect='auto',
           cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
plt.plot((msa_arr != "-").sum(0), color='black')
plt.xlim(-0.5,msa_arr.shape[1]-0.5)
plt.ylim(-0.5,msa_arr.shape[0]-0.5)
plt.colorbar(label="Sequence identity to query",)
plt.xlabel("Positions")
plt.ylabel("Sequences")


plt.subplot(1,2,2); plt.title("Predicted lDDT per position")
num_models = 0
for pkl_file in pkls:
	num_models += 1
	model_name = re.search(r'model_\d', pkl_file).group()
	with open(pkl_file, 'rb') as mpkl:
		resp = pickle.load(mpkl)
		plt.plot(resp['plddt'],label=model_name)

if homooligomer > 0:
  for n in range(homooligomer+1):
    x = n*(len(query_sequence)-1)
    plt.plot([x,x],[0,100],color="black")
plt.legend()
plt.ylim(0,100)
plt.ylabel("Predicted lDDT")
plt.xlabel("Positions")
plt.savefig("coverage_lDDT.pdf")
plt.savefig("coverage_lDDT.png")
#plt.show()
#
# plt.figure(figsize=(3*num_models,2), dpi=100)
# n = 0
#
# for pkl_file in pkls:
# 	model_name = re.search(r'model_\d', pkl_file).group()
# 	with open(pkl_file, 'rb') as mpkl:
# 		resp = pickle.load(mpkl)
# 		plt.subplot(1,num_models,n+1)
# 		plt.title(model_name)
# 		plt.imshow(resp["predicted_aligned_error"],label=model_name,cmap="bwr",vmin=0,vmax=30)
# 		plt.colorbar()
# 		n += 1
# plt.savefig("PAE.png")
#plt.show()
