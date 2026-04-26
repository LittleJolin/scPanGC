import scanpy as sc
import anndata as ad
import metacells as mc

def run_metacell_pipeline(adata, output_path):
    """
    原始 AnnData 数据，在同源组织背景下执行 MetaCell 降噪聚合。
    """
    print("\nStarting MetaCell construction...")
    celltypes = adata.obs['Celltype_new'].dropna().unique()
    diseases = adata.obs['Disease2'].dropna().unique()
    tissues = adata.obs['Tissue1'].dropna().unique()
    
    metacell_obj_list = []
    
    for celltype in celltypes:
        adata_tmp1 = adata[adata.obs['Celltype_new'] == celltype]
        for disease in diseases:
            adata_tmp2 = adata_tmp1[adata_tmp1.obs['Disease2'] == disease]
            valid_tissues = [t for t in tissues if t in adata_tmp2.obs['Tissue1'].values]
            
            for tissue in valid_tissues:
                adata_tmp3 = adata_tmp2[adata_tmp2.obs['Tissue1'] == tissue].copy()
                names_tmp = f"{celltype}-{disease}-{tissue}"
                
                if adata_tmp3.n_obs <= 50:
                    print(f"{names_tmp} NOT metacelled (n_obs={adata_tmp3.n_obs} <= 50)")
                    continue
                
                print(f"Processing MetaCell for: {names_tmp} (n_obs={adata_tmp3.n_obs})")
                mc.pl.divide_and_conquer_pipeline(adata_tmp3, target_metacell_size=50, random_seed=123456)
                metacells = mc.pl.collect_metacells(adata_tmp3, name=names_tmp+'.metacells')
                
                try:
                    mc.pl.compute_umap_by_features(metacells, max_top_feature_genes=1000, min_dist=2.0, random_seed=123456)
                except Exception as e:
                    print(f"  -> UMAP computation skipped for {names_tmp}: {e}")
                
                if '__name__' in metacells.uns:
                    del metacells.uns['__name__']
                    
                # 补全元数据标签
                metacells.obs['Celltype_new'] = celltype
                metacells.obs['Disease2'] = disease
                metacells.obs['Tissue1'] = tissue
                metacells.obs['Celltype3'] = celltype 
                metacells.obs['Disease'] = disease
                metacells.obs['Tissue'] = tissue
                metacell_obj_list.append(metacells)
                
    if len(metacell_obj_list) == 0:
        raise ValueError("No metacells were generated.")
        
    print(f"Concatenating {len(metacell_obj_list)} metacell chunks into a global atlas...")
    adata_metacell_global = ad.concat(metacell_obj_list, join="inner")
    adata_metacell_global.write_h5ad(output_path)
    
    return adata_metacell_global