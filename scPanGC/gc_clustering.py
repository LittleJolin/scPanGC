import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.stats as stats
import warnings

# 忽略计算过程中的除零等常规数学警告
warnings.filterwarnings('ignore')

def compute_gc_modules(adata_mc, deg_list, output_path, distance_threshold=70.0):
    """
    计算 Log2FC 矩阵并利用 Ward's 方法生成 Gene Clusters (GC)。
    """
    print("\nStarting GC Module computation...")
    
    # ========================================================
    # 🚀 1. 索引唯一化防线
    # 强制唯一化索引，彻底消除 "Observation names are not unique" 警告
    # ========================================================
    adata_mc.obs_names_make_unique()
    
    # ========================================================
    # 🚀 2. 分组与特征提取
    # ========================================================
    # 创建复合标签: 细胞类型 - 疾病 - 组织
    adata_mc.obs['C-D-T'] = adata_mc.obs['Celltype3'].astype(str) + '-' + \
                            adata_mc.obs['Disease'].astype(str) + '-' + \
                            adata_mc.obs['Tissue'].astype(str)
    
    # 取差异基因与元细胞矩阵基因的交集
    overlap_genes = [gene for gene in deg_list if gene in adata_mc.var_names]
    adata_sub = adata_mc[:, overlap_genes].copy()
    
    # 提取表达谱并按复合标签分组求均值
    obs_df = sc.get.obs_df(adata_sub, keys=overlap_genes)
    obs_df['C-D-T'] = adata_sub.obs['C-D-T'].values
    avg_score = obs_df.groupby('C-D-T')[overlap_genes].mean()
    
    # 伪计数平滑处理，防止取对数时遭遇 0
    avg_score = avg_score + 1.0 
    
    # ========================================================
    # 🚀 3. 计算 Log2FC 基准偏移矩阵
    # ========================================================
    cdt_list = avg_score.index.tolist()
    df_logfc = pd.DataFrame(index=cdt_list, columns=overlap_genes, dtype='float32')
    
    for cdt in cdt_list:
        parts = str(cdt).split('-')
        celltype, disease, tissue = parts[0], parts[1], parts[2]
        hc_name = f"{celltype}-HC-{tissue}"
        
        # 严格计算疾病组相对于同源健康对照(HC)的表达偏移
        if hc_name in avg_score.index:
            df_logfc.loc[cdt] = np.log2(avg_score.loc[cdt] / avg_score.loc[hc_name])
        else:
            df_logfc.loc[cdt] = 0.0 
            
    # 剔除全为 0 或 NaN 的无效行列
    df_logfc = df_logfc.dropna(axis=1, how='all').fillna(0)
    
    # ========================================================
    # 🚀 4. 层级聚类与模块切割
    # ========================================================
    # 对基因维度进行 Z-score 标准化
    z_matrix = stats.zscore(df_logfc.values, axis=0)
    
    # 采用 Ward's 方法构建聚类树
    distance_matrix = sch.distance.pdist(z_matrix.T, metric='euclidean')
    linkage_matrix = sch.linkage(distance_matrix, method='ward')
    
    # 按照设定的 distance_threshold (默认 70.0) 切割树，生成 GC
    cluster_labels = sch.fcluster(linkage_matrix, t=distance_threshold, criterion='distance')
    
    # 整理结果字典
    gc_result = pd.DataFrame({'Gene': df_logfc.columns, 'GC_Cluster': cluster_labels})
    
    # 自动建目录并保存结果
    out_dir = os.path.dirname(output_path)
    if out_dir: 
        os.makedirs(out_dir, exist_ok=True)
        
    gc_result.to_csv(output_path, index=False)
    print(f"Generated {len(set(cluster_labels))} Gene Clusters. Saved to {output_path}")
    
    return df_logfc, gc_result
