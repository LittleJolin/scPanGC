import os
import scanpy as sc
import pandas as pd
import warnings

def get_consensus_degs_from_raw(adata_raw, out_dir, min_cells=3, logfc_thresh=1.0, pval_thresh=0.05):
    """
    在原始单细胞矩阵中，严格限定“同源组织、同细胞亚群”背景，
    执行 Disease vs HC 的差异分析，并提取上调基因的并集。
    """
    print("Starting Differential Expression Analysis on RAW single-cell matrix...")
    
    # 获取唯一的分类标签
    celltypes = adata_raw.obs['Celltype_new'].dropna().unique()
    tissues = adata_raw.obs['Tissue1'].dropna().unique()
    diseases = [d for d in adata_raw.obs['Disease2'].dropna().unique() if d != 'HC']
    
    all_deg_list = []
    
    for celltype in celltypes:
        adata_ct = adata_raw[adata_raw.obs['Celltype_new'] == celltype]
        for tissue in tissues:
            adata_ct_tissue = adata_ct[adata_ct.obs['Tissue1'] == tissue]
            cells_hc = adata_ct_tissue[adata_ct_tissue.obs['Disease2'] == 'HC'].n_obs
            
            if cells_hc < min_cells:
                continue
            
            for disease in diseases:
                cells_disease = adata_ct_tissue[adata_ct_tissue.obs['Disease2'] == disease].n_obs
                if cells_disease >= min_cells:
                    adata_test = adata_ct_tissue[adata_ct_tissue.obs['Disease2'].isin([disease, 'HC'])].copy()
                    try:
                        sc.tl.rank_genes_groups(
                            adata_test, groupby='Disease2', groups=[disease], reference='HC', 
                            method='wilcoxon', use_raw=False
                        )
                        result = sc.get.rank_genes_groups_df(adata_test, group=disease)
                        sig_genes = result[(result['pvals_adj'] < pval_thresh) & (result['logfoldchanges'] >= logfc_thresh)]['names'].tolist()
                        if sig_genes:
                            all_deg_list.extend(sig_genes)
                            print(f"[{celltype} | {tissue} | {disease} vs HC] -> {len(sig_genes)} up-regulated DEGs.")
                    except Exception as e:
                        print(f"Skipped DEA for {celltype}-{tissue}-{disease} due to error: {e}")
                        
    final_deg_set = list(set(all_deg_list))
    # 剔除噪音基因
    final_deg_set = [g for g in final_deg_set if not str(g).startswith(('MT-', 'RPS', 'RPL'))]
    print(f"\nDEA completed! Total unique consensus up-regulated DEGs found: {len(final_deg_set)}")
    
    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame(final_deg_set, columns=['Gene']).to_csv(os.path.join(out_dir, 'consensus_DEGs.csv'), index=False)
    
    return final_deg_set