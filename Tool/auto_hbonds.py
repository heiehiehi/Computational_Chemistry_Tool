import os
import subprocess
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# 全局绘图风格设置 (学术风)
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2


def safe_decode(byte_data):
    if not byte_data:
        return ""
    try:
        return byte_data.decode('utf-8')
    except UnicodeDecodeError:
        # 如果 utf-8 解码失败，退回使用 gbk，并强行忽略无法识别的乱码字符
        return byte_data.decode('gbk', errors='ignore')


def win_to_wsl_path(win_path):
    """自动将 Windows 路径转换为 WSL 路径"""
    drive, tail = os.path.splitdrive(win_path)
    return f"/mnt/{drive[0].lower()}{tail.replace(os.sep, '/')}"


# ==========================================
# 步骤 0: 智能定位并加载全局配置
# ==========================================
# 智能寻路：获取当前脚本所在目录 (Tool) 的上一级目录 (工程根目录)
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

try:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
except FileNotFoundError:
    print(f"❌ 找不到 config.yaml 文件，请确保它存放在工程根目录: {PROJECT_ROOT}")
    exit(1)

BASE_DIR_WIN = Path(config['global']['base_dir_win'])
WSL_BASE_DIR = win_to_wsl_path(str(BASE_DIR_WIN))
TARGETS = " ".join(config['global']['targets'])
HIT_PREFIX = config['global']['hit_prefix']

# 聪明地借用 MMPBSA 配置中的环境和受体、配体组别
CONDA_ENV = config['mmpbsa']['conda_env']
GROUP_REC = str(config['mmpbsa']['group_receptor'])  # 默认 1
GROUP_LIG = str(config['mmpbsa']['group_ligand'])  # 默认 13

HBONDS_DIR_WIN = BASE_DIR_WIN / "HBonds"

# ==========================================
# 步骤 1: WSL 后台批量执行 GROMACS HBond
# ==========================================
wsl_bash_script_template = """
eval "$(conda shell.bash hook 2>/dev/null)" || source ~/miniconda3/etc/profile.d/conda.sh || source ~/anaconda3/etc/profile.d/conda.sh
conda activate __CONDA_ENV__

BASE_DIR="__WSL_BASE_DIR__"
HBONDS_DIR="${BASE_DIR}/HBonds"
TARGETS="__TARGETS__"
HIT_PREFIX="__HIT_PREFIX__"
GROUP_REC="__GROUP_REC__"
GROUP_LIG="__GROUP_LIG__"

echo "============================================================"
echo "🚀 步骤 1: 调用 WSL 批量执行 gmx hbond-legacy 计算氢键 (-tu ns)..."
echo "============================================================"

TOTAL_PROCESSED=0

for TARGET in $TARGETS; do
    if ls "${BASE_DIR}/${TARGET}"/${HIT_PREFIX}* 1> /dev/null 2>&1; then
        for HIT_DIR in "${BASE_DIR}/${TARGET}"/${HIT_PREFIX}*; do
            if [ -d "$HIT_DIR" ]; then
                HIT_NAME=$(basename "$HIT_DIR")

                # 建立对应的 靶点/Hit 文件夹
                OUT_DIR="${HBONDS_DIR}/${TARGET}/${HIT_NAME}"
                mkdir -p "$OUT_DIR"

                XVG_NUM="${OUT_DIR}/total_hbonds.xvg"
                XVG_LIFE="${OUT_DIR}/lifetime.xvg"

                # 断点续传机制
                if [ -f "$XVG_NUM" ]; then
                    echo "   ⏩ [已存在] 跳过 $TARGET -> $HIT_NAME"
                    continue
                fi

                echo "   ⚡ 正在计算 $TARGET -> $HIT_NAME 的氢键数量..."

                # 核心命令：换回 hbond-legacy，完美兼容 -life，同时加上 -tu ns 换算时间
                printf "%s\\n%s\\n" "$GROUP_REC" "$GROUP_LIG" | gmx hbond-legacy -f "${HIT_DIR}/md_fit.xtc" -s "${HIT_DIR}/md.tpr" -n "${HIT_DIR}/index.ndx" -num "$XVG_NUM" -life "$XVG_LIFE" -tu ns > /dev/null 2>&1

                if [ -f "$XVG_NUM" ]; then
                    echo "   ✅ 成功生成 $TARGET -> $HIT_NAME 的 XVG 文件！"
                    TOTAL_PROCESSED=$((TOTAL_PROCESSED + 1))
                else
                    echo "   ❌ 警告: $TARGET -> $HIT_NAME 未能生成 XVG！(请检查 index.ndx 是否存在组别)"
                fi
            fi
        done
    fi
done
echo "✅ WSL 计算部分执行完毕！(本次新生成 $TOTAL_PROCESSED 个体系)"
"""


def parse_hbond_xvg(file_path):
    """提取 XVG 文件中的有效数据，由于加了 -tu ns，第一列直接是时间(ns)，第二列是氢键数量"""
    data = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('@') or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    data.append((float(parts[0]), float(parts[1])))
    except Exception as e:
        print(f"读取 {file_path} 失败: {e}")
    # 注意这里列名改成了 ns
    return pd.DataFrame(data, columns=['Time (ns)', 'Hydrogen Bonds'])


def plot_hbond(df, target, hit_name, out_dir):
    """绘制高颜值学术氢键波动折线图"""
    fig, ax = plt.subplots(figsize=(12, 5), dpi=300)

    # 氢键绘图使用暖色调 (红色系)，与 RMSF 的冷色调区分开
    # 数据已经是 ns，直接使用 df['Time (ns)']，不用再除以 1000 了
    ax.plot(df['Time (ns)'], df['Hydrogen Bonds'], color='#D65F5F', linewidth=1.2, alpha=0.85)
    ax.fill_between(df['Time (ns)'], df['Hydrogen Bonds'], color='#D65F5F', alpha=0.25)

    ax.set_xlabel('Time (ns)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Number of Hydrogen Bonds', fontsize=14, fontweight='bold')
    ax.set_title(f'Protein-Ligand Hydrogen Bonds\n{target} - {hit_name}', fontsize=16, fontweight='bold')

    # 强制 Y 轴显示为整数（因为氢键数量不存在小数）
    ax.yaxis.get_major_locator().set_params(integer=True)

    ax.grid(axis='y', linestyle='--', alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()

    png_path = out_dir / f"{hit_name}_HBonds.png"
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    print("🚀 启动自动化 HBond (氢键) 提取与绘图脚本...\n")

    # 1. 组装并执行 WSL 脚本
    script = wsl_bash_script_template.replace("__CONDA_ENV__", CONDA_ENV) \
        .replace("__WSL_BASE_DIR__", WSL_BASE_DIR) \
        .replace("__TARGETS__", TARGETS) \
        .replace("__HIT_PREFIX__", HIT_PREFIX) \
        .replace("__GROUP_REC__", GROUP_REC) \
        .replace("__GROUP_LIG__", GROUP_LIG)
    try:
        process = subprocess.Popen(
            ['wsl', 'bash', '-i'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        )
        process.stdin.write(script.replace('\r\n', '\n').encode('utf-8'))
        process.stdin.close()
        for line in iter(process.stdout.readline, b''):
            print(safe_decode(line).rstrip())
        process.wait()
    except Exception as e:
        print(f"❌ WSL 运行异常: {e}")

    # 2. 回到 Windows 进行数据解析与绘图
    print("\n============================================================")
    print("🎨 步骤 2: 开始清洗数据，生成 CSV 与 高清氢键波动图...")
    print("============================================================")

    success_count = 0
    if HBONDS_DIR_WIN.exists():
        for target in config['global']['targets']:
            target_dir = HBONDS_DIR_WIN / target
            if not target_dir.exists(): continue

            for hit_dir in target_dir.iterdir():
                if hit_dir.is_dir() and hit_dir.name.startswith(HIT_PREFIX):
                    xvg_file = hit_dir / "total_hbonds.xvg"

                    if xvg_file.exists():
                        # 解析数据
                        df = parse_hbond_xvg(xvg_file)
                        if not df.empty:
                            # 导出 CSV
                            csv_path = hit_dir / f"{hit_dir.name}_HBonds.csv"
                            df.to_csv(csv_path, index=False)

                            # 绘制 PNG
                            plot_hbond(df, target, hit_dir.name, hit_dir)
                            print(f"   📊 [{target}] {hit_dir.name}: 图表与 CSV 转换完成！")
                            success_count += 1

    print("\n============================================================")
    print(f"🎉 全部任务完成！共处理了 {success_count} 个体系的 HBond 数据。")
    print(f"📂 产出结果已分发至: {HBONDS_DIR_WIN}")
    print("============================================================")