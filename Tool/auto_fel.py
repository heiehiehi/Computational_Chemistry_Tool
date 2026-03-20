import os
import subprocess
import yaml
import sys
import time
from pathlib import Path


# ==========================================
# 辅助函数
# ==========================================
def safe_decode(byte_data):
    if not byte_data: return ""
    try:
        return byte_data.decode('utf-8')
    except UnicodeDecodeError:
        return byte_data.decode('gbk', errors='ignore')


def win_to_wsl_path(win_path):
    drive, tail = os.path.splitdrive(win_path)
    return f"/mnt/{drive[0].lower()}{tail.replace(os.sep, '/')}"


# ==========================================
# 步骤 0: 加载配置
# ==========================================
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

try:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
except FileNotFoundError:
    print(f"❌ 找不到 config.yaml！请确保它在: {PROJECT_ROOT}")
    exit(1)

BASE_DIR_WIN = Path(config['global']['base_dir_win'])
WSL_BASE_DIR = win_to_wsl_path(str(BASE_DIR_WIN))

# 脚本路径修正：E:\GROMACS\Two Target\Script
SCRIPT_DIR_WIN = BASE_DIR_WIN / "Script"
WSL_SCRIPT_DIR = win_to_wsl_path(str(SCRIPT_DIR_WIN))

FEL_DIR_WIN = BASE_DIR_WIN / "FEL"
WSL_FEL_DIR = win_to_wsl_path(str(FEL_DIR_WIN))

TARGETS = config['global']['targets']
HIT_PREFIX = config['global']['hit_prefix']
CONDA_ENV = config['mmpbsa']['conda_env']


# ==========================================
# 核心函数：单体系处理逻辑 (WSL计算 + Windows画图)
# ==========================================
def process_single_hit(target, hit_name, hit_path_wsl):
    out_dir_win = FEL_DIR_WIN / target / hit_name
    out_dir_win.mkdir(parents=True, exist_ok=True)
    out_dir_wsl = win_to_wsl_path(str(out_dir_win))

    # 构建 Bash 命令 (1-9 步)
    bash_cmd = f"""
    eval "$(conda shell.bash hook 2>/dev/null)" || source ~/miniconda3/etc/profile.d/conda.sh || source ~/anaconda3/etc/profile.d/conda.sh
    conda activate {CONDA_ENV}
    cd "{out_dir_wsl}" || exit

    NDX="{hit_path_wsl}/index.ndx"
    TPR="{hit_path_wsl}/md.tpr"
    XTC="{hit_path_wsl}/md.xtc"

    # 1-3. PBC & Extraction
    printf "0\\n" | gmx trjconv -f "$XTC" -s "$TPR" -o complex.xtc -n "$NDX" -pbc mol -ur compact > /dev/null 2>&1
    printf "20\\n" | gmx trjconv -f complex.xtc -s "$TPR" -o complex2.xtc -n "$NDX" > /dev/null 2>&1
    printf "20\\n" | gmx trjconv -f complex.xtc -s "$TPR" -o complex2.gro -n "$NDX" -dump 0 > /dev/null 2>&1

    # 4. Fit
    printf "4\\n20\\n" | gmx trjconv -f complex2.xtc -s complex2.gro -o complex_fix.xtc -fit rot+trans -n "$NDX" > /dev/null 2>&1

    # 5. Covar (PCA)
    printf "4\\n4\\n" | gmx covar -f complex_fix.xtc -s complex2.gro -o eigenvalues.xvg -v eigenvectors.trr -xpma covapic.xpm -n "$NDX" > /dev/null 2>&1

    # 6. Proj
    printf "4\\n4\\n" | gmx anaeig -f complex_fix.xtc -s complex2.gro -v eigenvectors.trr -n "$NDX" -eig eigenvalues.xvg -first 1 -last 1 -proj pc1.xvg > /dev/null 2>&1
    printf "4\\n4\\n" | gmx anaeig -f complex_fix.xtc -s complex2.gro -v eigenvectors.trr -n "$NDX" -eig eigenvalues.xvg -first 2 -last 2 -proj pc2.xvg > /dev/null 2>&1

    # 7. Combine
    python3 "{WSL_SCRIPT_DIR}/pc_combine.py" pc1.xvg pc2.xvg output.xvg

    # 8. SHAM (300K)
    if [ -f "output.xvg" ]; then
        echo 1 | gmx sham -tsham 300 -nlevels 100 -f output.xvg -ls pc12_gibbs.xpm -g pc_12.log
    fi

    # 9. XPM2XYZ
    bash "{WSL_SCRIPT_DIR}/xpm2all.bsh" pc12_gibbs.xpm -xyz

    rm -f complex.xtc complex2.xtc
    """

    print(f"🚀 [WSL 计算中] {target} -> {hit_name}...")
    proc = subprocess.Popen(['wsl', 'bash'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, _ = proc.communicate(input=bash_cmd.encode('utf-8'))
    print(safe_decode(stdout))

    # --- 关键：计算完立刻启动 Windows 画图 ---
    xpm_file = out_dir_win / "pc12_gibbs.xpm"
    if xpm_file.exists():
        print(f"🎨 [Windows 绘图中] {target} -> {hit_name}...")
        draw_script = SCRIPT_DIR_WIN / "xpm2png.py"
        # 显式捕捉绘图脚本的错误输出
        res = subprocess.run([sys.executable, str(draw_script), "-ip", "yes", "-f", str(xpm_file)],
                             capture_output=True, text=True)
        if res.returncode == 0:
            print(f"✅ 图片已生成: {out_dir_win / 'pc12_gibbs.png'}")
        else:
            print(f"❌ 绘图失败内容: {res.stderr}")
    else:
        print(f"❌ 错误: WSL 未能生成 xpm 文件。")


# ==========================================
# 主程序入口
# ==========================================
if __name__ == "__main__":
    print("------------------------------------------------------------")
    print("🌟 PCA-FEL 增强版：即时流水线工作模式")
    print("------------------------------------------------------------")

    for target in TARGETS:
        target_dir_win = BASE_DIR_WIN / target
        if not target_dir_win.exists(): continue

        # 查找所有 Hit 子目录
        hits = [d for d in target_dir_win.iterdir() if d.is_dir() and d.name.startswith(HIT_PREFIX)]

        for hit in hits:
            # 断点续传检查
            if (FEL_DIR_WIN / target / hit.name / "pc12_gibbs.png").exists():
                print(f"⏩ [跳过] {target} -> {hit.name} (图片已存在)")
                continue

            process_single_hit(target, hit.name, win_to_wsl_path(str(hit)))

    print("\n🎉 全部流水线任务执行完毕！")