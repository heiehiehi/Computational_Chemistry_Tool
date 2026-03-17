import warnings

warnings.filterwarnings('ignore')

import os
import glob
import MDAnalysis as mda
import prolif as plf
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# 🌟 全局设置区：你只需要修改这里！
# ==========================================
# 输入你要扫描的靶点主目录（程序会自动提取最后面的 "JAK2" 或 "STAT3" 作为靶点名）
# 注意：Windows 路径前面最好加个 r，防止转义字符报错
BASE_DIR = r"E:\GROMACS\Two Target\STAT3"


# ==========================================

def process_single_hit(hit_path, out_dir, target_name, hit_name):
    """处理单个 Hit 文件夹的核心函数"""
    print(f"\n[{hit_name}] --------------------------------------------------")
    print(f"[{hit_name}] 开始处理...")

    # 1. 自动寻找需要的文件
    tpr_file = os.path.join(hit_path, "md.tpr")
    xtc_file = os.path.join(hit_path, "md_center.xtc")

    # 使用 glob 模糊匹配 pdb 文件 (只要是以 frame 开头的 pdb 都可以)
    pdb_files = glob.glob(os.path.join(hit_path, "frame*.pdb"))

    if not os.path.exists(tpr_file) or not os.path.exists(xtc_file):
        print(f"[{hit_name}] ❌ 警告：未找到 md.tpr 或 md_center.xtc，跳过此文件夹！")
        return
    if not pdb_files:
        print(f"[{hit_name}] ❌ 警告：未找到 frame_*.pdb 文件，跳过此文件夹！")
        return

    pdb_file = pdb_files[0]  # 取找到的第一个 pdb 文件

    try:
        # ==========================================
        # 步骤一：提取化学信息
        # ==========================================
        u_tpr = mda.Universe(tpr_file, xtc_file)
        lig = u_tpr.atoms.select_atoms("resname MOL")
        prot = u_tpr.atoms.select_atoms("protein")

        fp = plf.Fingerprint()
        fp.run(u_tpr.trajectory, lig, prot)  # 关闭单独的进度条，保持终端清爽

        df_raw = fp.to_dataframe()
        if df_raw.empty:
            print(f"[{hit_name}] ⚠️ 警告：该轨迹中未检测到任何有效相互作用！")
            return

        df_clean = df_raw.droplevel('ligand', axis=1)
        occupancy_series = df_clean.mean()

        # ==========================================
        # 步骤二：自动化纠正序号（利用 PDB 的真实序号）
        # ==========================================
        u_pdb = mda.Universe(pdb_file)
        prot_tpr = u_tpr.atoms.select_atoms("protein")
        prot_pdb = u_pdb.atoms.select_atoms("protein")

        auto_mapping = {}
        for r_tpr, r_pdb in zip(prot_tpr.residues, prot_pdb.residues):
            old_key = f"{r_tpr.resname}{r_tpr.resid}"
            new_val = f"{r_pdb.resname}{r_pdb.resid}"
            auto_mapping[old_key] = new_val

        # ==========================================
        # 步骤三：数据清洗与应用字典
        # ==========================================
        occ_df = occupancy_series.reset_index()
        occ_df.columns = ['Residue', 'Interaction', 'Occupancy']
        occ_df['Residue'] = occ_df['Residue'].astype(str).str.split('.').str[0]
        occ_df['Residue'] = occ_df['Residue'].replace(auto_mapping)

        pivot_df = occ_df.pivot_table(index='Residue', columns='Interaction', values='Occupancy', aggfunc='sum').fillna(
            0)
        pivot_df['Total_Occupancy'] = pivot_df.sum(axis=1)
        pivot_df = pivot_df[pivot_df['Total_Occupancy'] > 0.15]

        if pivot_df.empty:
            print(f"[{hit_name}] ⚠️ 警告：清洗后没有总占据率 > 0.15 的氨基酸。")
            return

        pivot_df = pivot_df.sort_values(by='Total_Occupancy', ascending=False).head(10)
        pivot_df = pivot_df.drop(columns=['Total_Occupancy'])

        # ==========================================
        # 步骤四：导出 Origin CSV 文件
        # ==========================================
        csv_path = os.path.join(out_dir, f"{hit_name}_Occupancy_Origin.csv")
        pivot_df.to_csv(csv_path)

        # ==========================================
        # 步骤五：绘制并保存图表
        # ==========================================
        interaction_colors_dict = {
            'VdWContact': '#4E79A7', 'Hydrophobic': '#59A14F', 'HBAcceptor': '#E15759',
            'HBDonor': '#F28E2B', 'PiStacking': '#B07AA1', 'CationPi': '#EDC948',
            'SaltBridge': '#76B7B2', 'MetalAcceptor': '#FF9DA7', 'PiCation': '#9C755F'
        }
        colors_to_use = [interaction_colors_dict.get(col, '#BAB0AC') for col in pivot_df.columns]

        fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
        pivot_df.plot(kind='bar', stacked=True, ax=ax, color=colors_to_use, width=0.75, edgecolor='black',
                      linewidth=0.5)

        ax.set_ylabel('Interaction Occupancy (Frequency)', fontsize=14, fontweight='bold')
        ax.set_xlabel('Key Amino Acid Residue', fontsize=14, fontweight='bold')
        ax.set_title(f'{target_name} - {hit_name} Protein-Ligand Interactions', fontsize=16, fontweight='bold')

        plt.xticks(rotation=45, ha='right', fontsize=12, fontweight='bold')
        plt.yticks(fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(title='Interaction Type', title_fontsize='13', fontsize='11',
                  bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, edgecolor='black')
        plt.tight_layout()

        png_path = os.path.join(out_dir, f"{hit_name}_Occupancy_Stacked.png")
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.close(fig)  # 极其重要：关闭图表释放内存，防止批处理时内存撑爆！

        print(f"[{hit_name}] ✅ 处理成功！结果已保存至: {out_dir}")

    except Exception as e:
        print(f"[{hit_name}] ❌ 处理失败，发生错误: {e}")


if __name__ == '__main__':
    # 获取当前执行程序的根目录
    CURRENT_DIR = os.getcwd()

    # 提取靶点名称 (比如把 "E:\GROMACS\Two Target\JAK2" 变成 "JAK2")
    TARGET_NAME = os.path.basename(os.path.normpath(BASE_DIR))

    print(f"🚀 启动全自动批处理任务...")
    print(f"📁 目标靶点: {TARGET_NAME}")
    print(f"🔍 扫描目录: {BASE_DIR}")

    # 检查基础目录是否存在
    if not os.path.exists(BASE_DIR):
        print(f"❌ 找不到目录: {BASE_DIR}，请检查路径是否正确！")
        exit()

    # 遍历该目录下的所有文件夹
    success_count = 0
    for folder_name in os.listdir(BASE_DIR):
        # 如果是以 "Hit" 开头，并且确实是个文件夹
        hit_path = os.path.join(BASE_DIR, folder_name)
        if folder_name.startswith("Hit") and os.path.isdir(hit_path):
            # 在当前程序目录下创建： 当前目录/JAK2/Hit11/
            out_hit_dir = os.path.join(CURRENT_DIR, TARGET_NAME, folder_name)
            os.makedirs(out_hit_dir, exist_ok=True)

            # 召唤核心函数去处理这个 Hit
            process_single_hit(hit_path, out_hit_dir, TARGET_NAME, folder_name)
            success_count += 1

    print("\n==================================================")
    print(f"🎉 批处理任务全部完成！共扫描了 {success_count} 个 Hit 文件夹。")
    print(f"📂 所有结果保存在: {os.path.join(CURRENT_DIR, TARGET_NAME)}")
    print("==================================================")