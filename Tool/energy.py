import subprocess
import sys


def safe_decode(byte_data):
    if not byte_data:
        return ""
    try:
        return byte_data.decode('utf-8')
    except UnicodeDecodeError:
        return byte_data.decode('gbk', errors='replace')


wsl_bash_script = """
# 加载 Conda 环境
eval "$(conda shell.bash hook 2>/dev/null)" || source ~/miniconda3/etc/profile.d/conda.sh || source ~/anaconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA 

BASE_DIR="/mnt/e/GROMACS/Two Target"
MMPBSA_IN="${BASE_DIR}/energy/mmpbsa.in"
OUT_BASE="${BASE_DIR}/MMPBSACal"

echo "============================================================"
echo "🔍 步骤 1: 检查环境与配置 (断点续传模式 🛡️)"
echo "============================================================"
if ! command -v gmx_MMPBSA &> /dev/null; then
    echo "❌ 致命错误: 找不到 gmx_MMPBSA 命令！"
    exit 1
fi

if ! command -v mpirun &> /dev/null; then
    echo "⚠️ 致命错误: 找不到 mpirun 命令！请确认你的环境中安装了 MPI。"
    exit 1
fi

if [ ! -f "$MMPBSA_IN" ]; then
    echo "❌ 致命错误: 找不到文件 $MMPBSA_IN"
    exit 1
fi

echo "📊 正在盘点计算任务..."
TOTAL_JOBS=0
for TARGET in STAT3 JAK2; do
    if ls "${BASE_DIR}/${TARGET}"/Hit* 1> /dev/null 2>&1; then
        for HIT_DIR in "${BASE_DIR}/${TARGET}"/Hit*; do
            if [ -d "$HIT_DIR" ]; then
                TOTAL_JOBS=$((TOTAL_JOBS + 1))
            fi
        done
    fi
done

echo "✅ 盘点完毕！总共有 $TOTAL_JOBS 个体系需要计算。"
echo ""
echo "▶️ 开始 【8核 MPI 并行 + 断点续传】 批量计算..."
echo ""

CURRENT_JOB=0
for TARGET in STAT3 JAK2; do
    if ls "${BASE_DIR}/${TARGET}"/Hit* 1> /dev/null 2>&1; then
        for HIT_DIR in "${BASE_DIR}/${TARGET}"/Hit*; do

            if [ -d "$HIT_DIR" ]; then
                CURRENT_JOB=$((CURRENT_JOB + 1))
                HIT_NAME=$(basename "$HIT_DIR")

                echo ""
                echo "============================================================"
                echo "⏳ [ 进度: $CURRENT_JOB / $TOTAL_JOBS ] 🎯 靶点: $TARGET | 体系: $HIT_NAME"
                echo "============================================================"

                OUT_DIR="${OUT_BASE}/${TARGET}/${HIT_NAME}"
                mkdir -p "$OUT_DIR"
                cd "$OUT_DIR" || exit

                # 🌟 【断点续传核心逻辑】
                if [ -f "FINAL_RESULTS_MMPBSA.dat" ]; then
                    echo "   ⏩ 发现 FINAL_RESULTS_MMPBSA.dat 文件！"
                    echo "   ✅ 此体系已计算完毕，跳过直接进行下一个..."
                    echo "------------------------------------------------------------"
                    continue
                fi

                echo "   📂 目录: $OUT_DIR"

                # 🌟 【拓扑文件智能修复逻辑】
                TOP_FILE="${HIT_DIR}/topol.top"
                if [ "$TARGET" == "JAK2" ]; then
                    echo "   🔧 检测到 JAK2 靶点，正在自动修复 PTR 缺失的二面角参数..."
                    # 【核心修复】：将 topol_fixed.top 生成在原始的 HIT_DIR 目录中，确保能找到包含的力场文件！
                    sed 's/\s\+4\s*$/    4   180.00   4.60240   2/' "${HIT_DIR}/topol.top" > "${HIT_DIR}/topol_fixed.top"
                    TOP_FILE="${HIT_DIR}/topol_fixed.top"
                    echo "   ✅ 拓扑文件修复完成，使用 ${HIT_DIR}/topol_fixed.top 进行计算"
                fi

                echo "   ⚡ 正在调用 8 个 CPU 物理核心狂暴计算中..."

                # 🌟 【8核提速 + 防报错 + 无弹窗 + 自动覆盖蓝屏残骸(-O)】
                mpirun --oversubscribe -np 8 gmx_MMPBSA -O -i "$MMPBSA_IN" \
                           -cs "${HIT_DIR}/md.tpr" \
                           -ct "${HIT_DIR}/md_fit.xtc" \
                           -ci "${HIT_DIR}/index.ndx" \
                           -cg 1 13 \
                           -cp "$TOP_FILE" \
                           -o FINAL_RESULTS_MMPBSA.dat \
                           -nogui

                echo "   ✅ [ $CURRENT_JOB / $TOTAL_JOBS ] $TARGET -> $HIT_NAME 计算完毕！"
                echo "------------------------------------------------------------"
            fi
        done
    fi
done

echo ""
echo "🎉🎉🎉 全部 $TOTAL_JOBS 个任务已圆满结束！ 🎉🎉🎉"
"""

if __name__ == "__main__":
    print("🚀 启动终极挂机版脚本 (断点续传 + 8核加速 + JAK2智能修复)...\n")
    try:
        process = subprocess.Popen(
            ['wsl', 'bash', '-i'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        )
        process.stdin.write(wsl_bash_script.encode('utf-8'))
        process.stdin.close()
        for line in iter(process.stdout.readline, b''):
            print(safe_decode(line).rstrip())
        process.wait()
    except Exception as e:
        print(f"❌ 运行异常: {e}")