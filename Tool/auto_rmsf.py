import warnings

warnings.filterwarnings('ignore')
import os
import glob
import yaml
import MDAnalysis as mda
import prolif as plf
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# 加载配置
with open("../config.yaml", "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

BASE_DIR_WIN = config['global']['base_dir_win']
TARGETS = config['global']['targets']
HIT_PREFIX = config['global']['hit_prefix']
LIG_RESNAME = config['prolif']['ligand_resname']
THRESHOLD = config['prolif']['occupancy_threshold']


def process_single_hit(hit_path, out_dir, target_name, hit_name):
    print(f"   [{hit_name}] 开始处理...")
    tpr_file = os.path.join(hit_path, "md.tpr")
    xtc_file = os.path.join(hit_path, "md_center.xtc")
    pdb_files = glob.glob(os.path.join(hit_path, "frame*.pdb"))

    if not os.path.exists(tpr_file) or not os.path.exists(xtc_file) or not pdb_files:
        print(f"   [{hit_name}] ❌ 缺少拓扑、轨迹或 PDB 文件，跳过！")
        return

    try:
        u_tpr = mda.Universe(tpr_file, xtc_file)
        lig = u_tpr.atoms.select_atoms(f"resname {LIG_RESNAME}")
        prot = u_tpr.atoms.select_atoms("protein")

        fp = plf.Fingerprint()
        fp.run(u_tpr.trajectory, lig, prot)
        df_raw = fp.to_dataframe()
        if df_raw.empty: return

        df_clean = df_raw.droplevel('ligand', axis=1)
        occupancy_series = df_clean.mean()

        u_pdb = mda.Universe(pdb_files[0])
        auto_mapping = {f"{r_t.resname}{r_t.resid}": f"{r_p.resname}{r_p.resid}"
                        for r_t, r_p in zip(prot.residues, u_pdb.atoms.select_atoms("protein").residues)}

        occ_df = occupancy_series.reset_index()
        occ_df.columns = ['Residue', 'Interaction', 'Occupancy']
        occ_df['Residue'] = occ_df['Residue'].astype(str).str.split('.').str[0].replace(auto_mapping)

        pivot_df = occ_df.pivot_table(index='Residue', columns='Interaction', values='Occupancy', aggfunc='sum').fillna(
            0)
        pivot_df['Total_Occupancy'] = pivot_df.sum(axis=1)
        pivot_df = pivot_df[pivot_df['Total_Occupancy'] > THRESHOLD].sort_values(by='Total_Occupancy',
                                                                                 ascending=False).head(10).drop(
            columns=['Total_Occupancy'])

        if pivot_df.empty: return

        # 输出 CSV
        pivot_df.to_csv(os.path.join(out_dir, f"{hit_name}_Occupancy_Origin.csv"))

        # 绘图
        colors_dict = {'VdWContact': '#4E79A7', 'Hydrophobic': '#59A14F', 'HBAcceptor': '#E15759', 'HBDonor': '#F28E2B'}
        colors = [colors_dict.get(c, '#BAB0AC') for c in pivot_df.columns]

        fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
        pivot_df.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.75, edgecolor='black', linewidth=0.5)
        ax.set_ylabel('Interaction Occupancy', fontsize=14, fontweight='bold')
        ax.set_title(f'{target_name} - {hit_name} Interactions', fontsize=16, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{hit_name}_Occupancy_Stacked.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)

        print(f"   ✅ [{hit_name}] 分析完成！")
    except Exception as e:
        print(f"   ❌ [{hit_name}] 错误: {e}")


if __name__ == '__main__':
    print(f"🚀 启动 Auto-ProLIF 批处理任务...\n")
    CURRENT_DIR = os.getcwd()

    for target in TARGETS:
        target_path = os.path.join(BASE_DIR_WIN, target)
        if not os.path.exists(target_path): continue

        print(f"📁 正在处理靶点: {target}")
        for folder_name in os.listdir(target_path):
            hit_path = os.path.join(target_path, folder_name)
            if folder_name.startswith(HIT_PREFIX) and os.path.isdir(hit_path):
                out_hit_dir = os.path.join(CURRENT_DIR, target, folder_name)
                os.makedirs(out_hit_dir, exist_ok=True)
                process_single_hit(hit_path, out_hit_dir, target, folder_name)