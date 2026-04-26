import scanpy as sc
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.stats as stats

def compute_gc_modules(adata_mc, deg_list, output_path, distance_threshold=70.0):
    """
    计算 Log2FC 矩阵并利用 Ward's 方法生成 Gene Clusters (GC)。
    """
    print("\nStarting GC Module computation...")
    # 创建复合标签
    adata_mc.obs['C-D-T'] = adata_mc.obs['Celltype3'].astype(str) + '-' + \
                            adata_mc.obs['Disease'].astype(str) + '-' + \
                            adata_mc.obs['Tissue'].astype(str)
    
    overlap_genes = [gene for gene in deg_list if gene in adata_mc.var_names]
    adata_sub = adata_mc[:, overlap_genes].copy()
    
    obs_df = sc.get.obs_df(adata_sub, keys=overlap_genes)
    obs_df['C-D-T'] = adata_sub.obs['C-D-T'].values
    avg_score = obs_df.groupby('C-D-T')[overlap_genes].mean()
    avg_score = avg_score + 1.0 # 伪计数处理
    
    cdt_list = avg_score.index.tolist()
    df_logfc = pd.DataFrame(index=cdt_list, columns=overlap_genes, dtype='float32')
    
    for cdt in cdt_list:
        parts = str(cdt).split('-')
        celltype, disease, tissue = parts[0], parts[1], parts[2]
        hc_name = f"{celltype}-HC-{tissue}"
        
        if hc_name in avg_score.index:
            df_logfc.loc[cdt] = np.log2(avg_score.loc[cdt] / avg_score.loc[hc_name])
        else:
            df_logfc.loc[cdt] = 0.0 
            
    df_logfc = df_logfc.dropna(axis=1, how='all').fillna(0)
    
    # 执行标准化与层级聚类
    z_matrix = stats.zscore(df_logfc.values, axis=0)
    distance_matrix = sch.distance.pdist(z_matrix.T, metric='euclidean')
    linkage_matrix = sch.linkage(distance_matrix, method='ward')
    cluster_labels = sch.fcluster(linkage_matrix, t=distance_threshold, criterion='distance')
    
    gc_result = pd.DataFrame({'Gene': df_logfc.columns, 'GC_Cluster': cluster_labels})
    gc_result.to_csv(output_path, index=False)
    
    return df_logfc, gc_result