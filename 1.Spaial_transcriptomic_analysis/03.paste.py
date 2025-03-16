import math
import time
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import style
import paste as past


input_path="cs9_embryo_annotation.h5ad"
sample= pd.read_csv("sample_name.csv")

##Read data and create AnnData object
cs9=sc.read_h5ad(input_path)
cs9.obs['barcode']=cs9.obs.index.tolist()
cs9.X=cs9.raw.X
sc.pp.calculate_qc_metrics(cs9, inplace=True)
cs9.obsm['spatial'] = np.array(cs9.obs[['newx1','newy1']])
slice_names=sample['V1'].tolist()

def load_slices(h5ad,slice_names):
    slices = []  
    for slice_name in slice_names:
        slice_i = h5ad[h5ad.obs['sample'] == slice_name,:]
        # Preprocess slices
        #sc.pp.filter_genes(slice_i, min_counts = 15)
        #sc.pp.filter_cells(slice_i, min_counts = 100)
        slices.append(slice_i)
    return slices

slices = load_slices(cs9,slice_names)

##Pairwise Alignment
pis=[]
for i in range(len(slices)-1):
    pi_temp = past.pairwise_align(slices[i], slices[i+1])
    pis.append(pi_temp)

##Sequential pairwise slice alignment plots
new_slices = past.stack_slices_pairwise(slices, pis)


## plot the slices in 3-D.
# scale the distance between layers
z_scale = 0.5

TrD_module= pd.DataFrame(0, columns=['x','y','z','barcode'],index=range(0))
for i,L in enumerate(new_slices):
    values = []
    for x,y in L.obsm['spatial']:
        values.append([x, y, i*z_scale])
        temp= pd.DataFrame(values, columns=['x','y','z'])
    temp['barcode']=L.obs['barcode'].tolist()
    TrD_module=pd.concat([TrD_module,temp],axis=0)

TrD_module.to_csv(${output_path},index=True)