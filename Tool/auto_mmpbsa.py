import os
import subprocess
import yaml
from pathlib import Path


def safe_decode(byte_data):
    return byte_data.decode('utf-8') if byte_data else ""


def win_to_wsl_path(win_path):
    """自动将 Windows 路径 (E:\xxx) 转换为 WSL 路径 (/mnt/e/xxx)"""
    drive, tail = os.path.splitdrive(win_path)
    return f"/mnt/{drive[0].lower()}{tail.replace(os.sep, '/')}"


# 加载配置文件
with open("../config.yaml", "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

# 提取配置参数
WIN_BASE_DIR = config['global']['base_dir_win']
WSL_BASE_DIR = win_to_wsl_path(WIN_BASE_DIR)
TARGETS = " ".join(config['global']['targets'])
HIT_PREFIX = config['global']['hit_prefix']

CONDA_ENV = config['mmpbsa']['conda_env']
MMPBSA_IN_WSL = f"{WSL_BASE_DIR}/{config['mmpbsa']['mmpbsa_in_rel']}"
CORES = config['mmpbsa']['mpi_cores']
CG_REC = config['mmpbsa']['group_receptor']
CG_LIG = config['mmpbsa']['group_ligand']
FIX_TARGETS = config['mmpbsa'].get('auto_fix_targets', [])

wsl_bash_script_template = """
eval "$(conda shell.bash hook 2>/dev/null)" || source ~/miniconda3/etc/profile.d/conda.sh || source ~/anaconda3/etc/profile.d/conda.sh
conda activate __CONDA_ENV__

BASE_DIR="__WSL_BASE_DIR__"
MMPBSA_IN="__MMPBSA_IN__"
OUT_BASE="${BASE_DIR}/MMPBSACal"
TARGETS="__TARGETS__"
HIT_PREFIX="__HIT_PREFIX__"

echo "============================================================"
echo "🚀 启动 MMPBSA 自动化流水线 (配置已加载)"
echo "============================================================"

TOTAL_JOBS=0
for TARGET in $TARGETS; do
    if ls "${BASE_DIR}/${TARGET}"/${HIT_PREFIX}* 1> /dev/null 2>&1; then
        for HIT_DIR in "${BASE_DIR}/${TARGET}"/${HIT_PREFIX}*; do
            if [ -d "$HIT_DIR" ]; then TOTAL_JOBS=$((TOTAL_JOBS + 1)); fi
        done
    fi
done

echo "📊 发现总任务数: $TOTAL_JOBS"

CURRENT_JOB=0
for TARGET in $TARGETS; do
    if ls "${BASE_DIR}/${TARGET}"/${HIT_PREFIX}* 1> /dev/null 2>&1; then
        for HIT_DIR in "${BASE_DIR}/${TARGET}"/${HIT_PREFIX}*; do
            if [ -d "$HIT_DIR" ]; then
                CURRENT_JOB=$((CURRENT_JOB + 1))
                HIT_NAME=$(basename "$HIT_DIR")
                OUT_DIR="${OUT_BASE}/${TARGET}/${HIT_NAME}"
                mkdir -p "$OUT_DIR"
                cd "$OUT_DIR" || exit

                if [ -f "FINAL_RESULTS_MMPBSA.dat" ]; then
                    echo "   ⏩ [已完成] 跳过 $TARGET -> $HIT_NAME"
                    continue
                fi

                echo "------------------------------------------------------------"
                echo "⏳ [$CURRENT_JOB / $TOTAL_JOBS] 计算体系: $TARGET | $HIT_NAME"

                TOP_FILE="${HIT_DIR}/topol.top"

                # 动态注入智能修复逻辑 (根据 YAML 配置)
                __FIX_LOGIC__

                echo "   ⚡ 启动 __CORES__ 核狂暴计算..."
                mpirun --oversubscribe -np __CORES__ gmx_MMPBSA -O -i "$MMPBSA_IN" \
                           -cs "${HIT_DIR}/md.tpr" -ct "${HIT_DIR}/md_fit.xtc" -ci "${HIT_DIR}/index.ndx" \
                           -cg __CG_REC__ __CG_LIG__ -cp "$TOP_FILE" -o FINAL_RESULTS_MMPBSA.dat -nogui < /dev/null

                if [ -f "FINAL_RESULTS_MMPBSA.dat" ]; then
                    echo "   ✅ 成功！"
                else
                    echo "   ❌ 失败！"
                fi
            fi
        done
    fi
done
echo "🎉 所有 MMPBSA 任务执行完毕！"
"""

# 生成 Bash 中需要修复 Topology 的条件语句
fix_logic_bash = ""
if FIX_TARGETS:
    targets_condition = " -o ".join([f'[ "$TARGET" == "{t}" ]' for t in FIX_TARGETS])
    fix_logic_bash = f"""
                if [ {targets_condition} ]; then
                    echo "   🔧 [智能修复] 修改拓扑中缺失的参数..."
                    sed 's/\\s\\+4\\s*$/    4   180.00   4.60240   2/' "${{HIT_DIR}}/topol.top" > "${{HIT_DIR}}/topol_fixed.top"
                    TOP_FILE="${{HIT_DIR}}/topol_fixed.top"
                fi
    """

# 替换占位符
script = wsl_bash_script_template.replace("__CONDA_ENV__", CONDA_ENV) \
    .replace("__WSL_BASE_DIR__", WSL_BASE_DIR) \
    .replace("__MMPBSA_IN__", MMPBSA_IN_WSL) \
    .replace("__TARGETS__", TARGETS) \
    .replace("__HIT_PREFIX__", HIT_PREFIX) \
    .replace("__CORES__", str(CORES)) \
    .replace("__CG_REC__", str(CG_REC)) \
    .replace("__CG_LIG__", str(CG_LIG)) \
    .replace("__FIX_LOGIC__", fix_logic_bash)

if __name__ == "__main__":
    try:
        process = subprocess.Popen(['wsl', 'bash', '-i'], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        process.stdin.write(script.replace('\r\n', '\n').encode('utf-8'))
        process.stdin.close()
        for line in iter(process.stdout.readline, b''):
            print(safe_decode(line).rstrip())
        process.wait()
    except Exception as e:
        print(f"❌ 运行异常: {e}")