import os
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

try:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
except FileNotFoundError:
    print(f"找不到配置文件，请确保 config.yaml 存放在工程根目录: {PROJECT_ROOT}")
    exit(1)

BASE_DIR_WIN = Path(config['global']['base_dir_win'])
TARGETS = config['global']['targets']
HIT_PREFIX = config['global']['hit_prefix']
LIGAND_RESNAME = config['prolif'].get('ligand_resname', 'MOL')

MMPBSA_DIR_WIN = BASE_DIR_WIN / "MMPBSACal"
DRESULT_DIR_WIN = BASE_DIR_WIN / "Dresult"


def format_residue(raw_res):
    parts = raw_res.split(':')
    if len(parts) >= 2:
        res_name = parts[-2]
        res_num = parts[-1]
        return f"{res_name}{res_num}"
    return raw_res


def parse_decomp_dat(file_path):
    data_lines = []
    in_deltas = False
    capture = False
    total_idx = -1

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()

                if line == "DELTAS:":
                    in_deltas = True
                    continue

                if in_deltas and line == "Total Energy Decomposition:":
                    capture = True
                    continue

                if capture:
                    if not line:
                        break

                    parts = [p.strip() for p in line.split(',')]

                    if "Residue" in parts and "TOTAL" in parts:
                        total_idx = parts.index("TOTAL")
                        continue

                    if parts[0] == "" and "Avg." in parts:
                        continue

                    if total_idx != -1 and len(parts) > total_idx + 2:
                        raw_res = parts[0]
                        try:
                            total_ene = float(parts[total_idx])
                            total_err = float(parts[total_idx + 2])
                            data_lines.append((raw_res, total_ene, total_err))
                        except ValueError:
                            pass
    except Exception as e:
        print(f"读取 {file_path} 失败: {e}")
        return pd.DataFrame()

    df = pd.DataFrame(data_lines, columns=['Raw_Residue', 'TOTAL_Energy', 'TOTAL_Error'])
    return df


def plot_top6(df, target, hit_name, out_dir):
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)

    ax.bar(df['Formatted_Residue'], df['TOTAL_Energy'], yerr=df['TOTAL_Error'],
           capsize=5, color='#4C72B0', edgecolor='black', linewidth=1)

    ax.set_xlabel('Key Amino Acid Residue', fontsize=12)
    ax.set_ylabel('Energy Contribution (kcal/mol)', fontsize=12)
    ax.set_title(f'Top 6 Binding Energy Contributions\n{target} - {hit_name}', fontsize=14)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linestyle='--', alpha=0.6)

    plt.xticks(rotation=45, ha='right', fontsize=11)
    plt.yticks(fontsize=11)
    plt.tight_layout()

    png_path = out_dir / f"{hit_name}_Top6_Energy.png"
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    print("启动结合能分解提取与绘图脚本")

    success_count = 0
    summary_data = []

    if MMPBSA_DIR_WIN.exists():
        for target in TARGETS:
            target_in_dir = MMPBSA_DIR_WIN / target
            target_out_dir = DRESULT_DIR_WIN / target

            if not target_in_dir.exists():
                continue

            for hit_dir in target_in_dir.iterdir():
                if hit_dir.is_dir() and hit_dir.name.startswith(HIT_PREFIX):
                    dat_file = hit_dir / "FINAL_DECOMP_MMPBSA.dat"

                    if dat_file.exists():
                        df = parse_decomp_dat(dat_file)

                        if not df.empty:
                            df = df[~df['Raw_Residue'].str.contains(LIGAND_RESNAME, case=False, na=False)]

                            df['Formatted_Residue'] = df['Raw_Residue'].apply(format_residue)

                            out_hit_dir = target_out_dir / hit_dir.name
                            out_hit_dir.mkdir(parents=True, exist_ok=True)

                            all_res_df = df[['Formatted_Residue', 'TOTAL_Energy', 'TOTAL_Error']]
                            all_csv_path = out_hit_dir / f"{hit_dir.name}_All_Residues_Decomp.csv"
                            all_res_df.to_csv(all_csv_path, index=False)

                            top6_df = df.nsmallest(6, 'TOTAL_Energy')
                            top6_csv_path = out_hit_dir / f"{hit_dir.name}_Top6_Residues.csv"
                            top6_df.to_csv(top6_csv_path, index=False)

                            plot_top6(top6_df, target, hit_dir.name, out_hit_dir)

                            top1_res = top6_df.iloc[0]
                            summary_data.append({
                                'Target': target,
                                'Hit': hit_dir.name,
                                'Strongest_Residue': top1_res['Formatted_Residue'],
                                'Energy': top1_res['TOTAL_Energy'],
                                'Error': top1_res['TOTAL_Error']
                            })

                            print(f"成功处理 {target} -> {hit_dir.name}")
                            success_count += 1

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        DRESULT_DIR_WIN.mkdir(parents=True, exist_ok=True)
        summary_csv_path = DRESULT_DIR_WIN / "Strongest_Residue_Summary.csv"
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"已生成全局最强氨基酸统计表: {summary_csv_path}")

    print(f"任务完成，共处理了 {success_count} 个体系的数据。")
    print(f"结果已保存至: {DRESULT_DIR_WIN}")