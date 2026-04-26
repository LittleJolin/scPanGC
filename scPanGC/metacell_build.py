import scanpy as sc
import anndata as ad
import metacells as mc

def run_metacell_pipeline(adata, output_path):
    """
    接收原始 AnnData 数据，在同源组织背景下执行 MetaCell 降噪聚合。
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
                # ========================================================
                # 🚀 新增：补齐 Metacells 强制要求的基因掩码 (Masks)
                # ========================================================
                for required_mask in ['lateral_gene', 'noisy_gene']:
                    if required_mask not in adata_tmp3.var.columns:
                        adata_tmp3.var[required_mask] = False
                # ========================================================

                # 继续执行后续 MetaCell 流程
                mc.pl.divide_and_conquer_pipeline(adata_tmp3, target_metacell_size=50, random_seed=123456)
                
                # 【此处已修复】: 补充了新版 metacells 强制要求的 random_seed 参数
                metacells = mc.pl.collect_metacells(adata_tmp3, name=names_tmp+'.metacells', random_seed=123456)
                
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
        

import os # 确保文件开头或者这里引入了 os 模块

    if len(metacell_obj_list) == 0:
        raise ValueError("No metacells were generated. Please check your data quality.")
        
    print(f"\nConcatenating {len(metacell_obj_list)} metacell chunks into a global atlas...")
    adata_metacell_global = ad.concat(metacell_obj_list, join="inner")
    
    # 🚀 修复1：强制让所有的元细胞拥有唯一的 barcode 名字 (避免串号)
    adata_metacell_global.obs_names_make_unique()
    
    # 🚀 修复2：自动检查并创建输出文件夹，防止 FileNotFoundError
    out_dir = os.path.dirname(output_path)
    if out_dir: # 如果路径包含文件夹层级
        os.makedirs(out_dir, exist_ok=True)
        
    # 安全保存
    adata_metacell_global.write_h5ad(output_path)
    print(f"Global Metacell construction completed. Saved to {output_path}")
    
    return adata_metacell_global
