import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt
import anndata as ad
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

sample=["M2-1","M2-4","M2-5","M2-6","M2-7","M2-8","M2-9","M2-10",
        "M1-1","M1-2","M1-3","M1-4","M1-5","M1-6","M1-7","M1-8","M1-9","M1-10",
         "M3-1","M3-2","M3-3","M3-4","M3-5","M3-6","M3-7","M3-8","M3-9","M3-10",
         "M5-1","M5-2","M5-3","M5-4","M5-5","M5-6","M5-7","M5-8","M5-9","M5-10"]
adata_rna_list=[]
adata_atac_list=[]

for x in range(len(sample)):
 target=sample[x]
 print(target) 
 atac_path = '/database/wangrong/Results/0712_ATAC+RNA/' + target + '/outs/filtered_feature_bc_matrix'
 atac_annot_path = '/database/wangrong/Results/0712_ATAC+RNA/' + target + '/outs/atac_peak_annotation.tsv'
 feature_linkage_path = '/database/wangrong/Results/0712_ATAC+RNA/' + target + '/outs/analysis/feature_linkage/feature_linkage.bedpe'
 loom_path = '/database/wangrong/Results/0712_ATAC+RNA/' + target + '/velocyto/' + target + '.loom'
 annot_path = '/home/huangdingli/sample_metadata.csv'


 adata_rna = scv.read(loom_path, cache=True)
 adata_rna.obs_names = [x.replace(':', '_')[:-1] + '-1' for x in adata_rna.obs_names]  # 加上'-1'后缀
 adata_rna.var_names_make_unique()
 sc.pp.filter_cells(adata_rna, min_counts=1000)
 sc.pp.filter_cells(adata_rna, max_counts=20000)

 cell_annot = pd.read_csv(annot_path, sep=',', index_col=0)
 cell_annot_1= cell_annot[cell_annot["orig.ident"] == target]
 shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, cell_annot_1.index))
 adata_rna = adata_rna[shared_cells, ]
 adata_rna
 #print(cell_annot.head())
 cell_annot_1.index
 shared_cells

 cell_annot_2 = cell_annot_1.loc[shared_cells]
 adata_rna.obs['celltype'] = cell_annot_2['celltype']
 adata_rna.obs['orig.ident'] = cell_annot_2['orig.ident']
 adata_rna.obs['Timepoints'] = cell_annot_2['Timepoints']
 adata_rna.obs['Participants'] = cell_annot_2['Participants']
 adata_rna.obs['annotation'] = cell_annot_2['annotation']
 adata_rna

 adata_atac = sc.read_10x_mtx(atac_path,
                             var_names='gene_symbols',
                             cache=True, gex_only=False)
 adata_atac = adata_atac[:, adata_atac.var['feature_types'] == "Peaks"]
 adata_atac = mv.aggregate_peaks_10x(adata_atac,
                                    atac_annot_path,
                                    feature_linkage_path,
                                    verbose=True)
 adata_atac.obs_names = [target + '_' + x  for x in adata_atac.obs_names] 
 plt.hist(adata_atac.X.sum(1), bins=100, range=(0, 40000))
 plt.show()
 sc.pp.filter_cells(adata_atac, min_counts=1500)
 sc.pp.filter_cells(adata_atac, max_counts=25000)

 adata_rna_list.append(adata_rna)
 adata_atac_list.append(adata_atac)
 
all_adata_rna=ad.concat(adata_rna_list,merge = "same")
all_adata_rna
all_adata_atac=ad.concat(adata_atac_list,merge = "same")
all_adata_atac

all_adata_rna.write("/database/wangrong/Results/0712_ATAC+RNA/python/10.28.all_adata_rna.h5ad")
all_adata_atac.write("/database/wangrong/Results/0712_ATAC+RNA/python/10.28.all_adata_atac.h5ad")