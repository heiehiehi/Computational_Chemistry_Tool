# Tool工具集合

## 该工具本作者使用在Win下的WSL2平台，系统为Ubuntu22，个人使用请依据自己数据进行修改（这些工具只是整合了一些流程）

# Auto-CADD-Pipeline 自动化分析

## 1. 简介
本项目是我在日常处理分子动力学（MD）模拟数据时，随手编写并整合的一套常用后处理分析工具脚本。

主要包含几个常用的分析模块：
1. **MM/PBSA 结合自由能批量计算**
2. **ProLIF 蛋白-配体相互作用指纹分析**
3. **RMSF 均方根涨落提取与绘图**
4. **氢键数目 提取与绘图**


不用去改任何 Python 源代码，只需要在全局配置文件中填好你的路径和靶点名，就能全自动跑批处理。

## 2. 运行环境要求
本项目基于 Windows + WSL2 (Linux) 混合运行架构，请确保你的电脑已具备以下基础环境：

* **Windows 端**：
  * Python 3.11+
  * 安装依赖库：`pip install pyyaml MDAnalysis prolif pandas matplotlib`
* **WSL 端 (如 Ubuntu 22.04)**：
  * GROMACS 5.x+ （确保 `gmx` 命令可用）
  * OpenMPI （用于多核加速，安装命令：`sudo apt install openmpi-bin`）
  * 安装并配置好包含 `gmx_MMPBSA` 及 AmberTools 的 Conda 虚拟环境。

## 3. 全局配置文件 (`config.yaml`)
这是整个工具箱的大脑。所有的路径、靶点名、组别编号都集中在这里管理。文件内部已经写好了详细的注释，**直接对照注释修改成你自己的信息即可**。

***

## 🔬 工具一：`auto_mmpbsa.py` (MM/PBSA 结合自由能批量计算)

### 1. 工具简介
这个脚本主要用来在后台（WSL 环境）全自动、批量地运行 `gmx_MMPBSA`，计算配体与受体的结合自由能。

**它的核心作用是帮你省去重复敲命令的麻烦，并内置了几个实用功能：**
* **多核加速**：自动调用 MPI 分配多个 CPU 核心同时跑，大幅缩短长轨迹的计算时间。
* **断点续传**：跑一半如果电脑死机、断电，或者你手动掐断，下次重新运行脚本时，它会自动跳过已经算完（已生成 `FINAL_RESULTS_MMPBSA.dat`）的体系，接着往下算。
* **拓扑自动修复**：针对像 `JAK2` 这种带有磷酸化残基（PTR）的体系，底层的 `ParmEd` 翻译拓扑时经常会因为找不到 `impropers` 或 `func 4` 参数报错。脚本会根据 `config.yaml` 的配置，自动对 `topol.top` 进行修改补全，直接规避报错。

### 2. 目录结构与依赖文件
要让脚本正常运行，你需要按照以下层级准备好原始文件。程序会自动去读取 `config.yaml` 中设置的 `base_dir_win` 目录：

```text
E:\GROMACS\Two Target\          <-- 你的工程根目录 (config.yaml 中的 base_dir_win)
├── energy\                     
│   └── mmpbsa.in               <-- 【必须】MMPBSA 的输入参数控制文件 (需指定好力场等)
├── STAT3\                      <-- 靶点名称 (config 中的 targets)
│   ├── Hit1\                   <-- 具体的复合物体系 (以 config 中的 hit_prefix 开头)
│   │   ├── md_fit.xtc          <-- 【必须】已经去除 PBC (周期性边界条件) 并居中叠合的轨迹文件
│   │   ├── md.tpr              <-- 【必须】GROMACS 的 MD 运行输入文件
│   │   ├── topol.top           <-- 【必须】完整的复合物拓扑文件
│   │   └── index.ndx           <-- 【必须】索引文件 (需包含受体和配体的独立组别)
│   └── Hit2\
└── MMPBSACal\                  <-- 【运行后自动生成】存放所有体系计算结果的总目录
```

### 3. 操作指南
1. **核对配置**：打开 `config.yaml`，确保 `mmpbsa` 节点下的 `conda_env`（环境名）、`mpi_cores`（CPU核心数）、`group_receptor` 和 `group_ligand`（受体配体的索引号）都填写正确。
2. **检查文件**：确认所有待计算的 Hit 目录中，那 4 个【必须】的基础文件都齐了。
3. **一键运行**：在 Windows 的终端（CMD、PowerShell 或 PyCharm 终端）中，直接执行：
   ```bash
   python auto_mmpbsa.py
   ```
4. **收取结果**：程序处于静默无弹窗模式，你可以直接挂机。跑完后，去工程根目录新生成的 `MMPBSACal` 文件夹里，按靶点和 Hit 名字找对应的 `FINAL_RESULTS_MMPBSA.dat` 结果文件即可。

***

## 📊 工具二：`auto_prolif.py` (ProLIF 蛋白-配体相互作用指纹分析)

### 1. 工具简介
这个脚本基于 `MDAnalysis` 和 `ProLIF` 库，用来批量提取 MD 轨迹中蛋白与配体之间的动态非共价相互作用（如氢键、盐桥、π-π 堆叠、疏水作用等）。

**它的核心作用是帮你把枯燥的轨迹转化为直观的图表，并解决了 GROMACS 序号错位的痛点：**
* **自动映射真实残基序号**：GROMACS 生成 `.tpr` 拓扑时会把残基序号打乱重排。脚本会读取你提供的一帧 `.pdb` 结构，自动将内部错位的序号完美映射回晶体学真实的残基编号，保证出图时的标签绝对准确。
* **数据提纯与过滤**：程序会计算整个 MD 过程中各相互作用的发生频率（占据率），并根据你在 `config.yaml` 中设置的阈值（如 `0.15` 代表 15%）自动过滤掉低频噪音，只保留贡献最大的 Top 10 关键氨基酸。
* **双轨输出**：每个体系跑完后，会直接生成一张带有标准配色、可用于文章的高清堆叠柱状图（`.png`），同时导出一份清洗好的干净数据表（`.csv`），方便你随时拖入 Origin 软件进行二次作图。

### 2. 目录结构与依赖文件
脚本运行时会读取 `config.yaml` 里的 `base_dir_win` 目录，并在**当前脚本所在的文件夹下**自动创建同名靶点目录来存放结果。

原始数据请确保按照以下结构存放：

```text
E:\GROMACS\Two Target\          <-- 你的工程根目录 (config.yaml 中的 base_dir_win)
├── STAT3\                      <-- 靶点名称 (config 中的 targets)
│   ├── Hit1\                   <-- 具体的复合物体系 (以 config 中的 hit_prefix 开头)
│   │   ├── md_center.xtc       <-- 【必须】已去除 PBC 并居中的轨迹文件
│   │   ├── md.tpr              <-- 【必须】MD 运行输入拓扑文件
│   │   └── frame_0.pdb         <-- 【必须】任意一帧的 PDB 文件 (脚本靠它来获取真实的残基序号)
│   └── Hit2\
└── 当前运行目录\                 <-- 【运行后自动生成】
    └── STAT3\                  
        ├── Hit1\
        │   ├── Hit1_Occupancy_Origin.csv    <-- 清洗好的数据源
        │   └── Hit1_Occupancy_Stacked.png   <-- 直出的堆叠柱状图
        └── Hit2\
```

### 3. 操作指南
1. **核对配置**：打开 `config.yaml`，检查 `prolif` 节点下的 `ligand_resname`（你的小分子在拓扑里的残基名，比如 `MOL`）和 `occupancy_threshold`（你想要的过滤阈值，默认 `0.15`）。
2. **检查文件**：确认每个待处理的 Hit 目录中都有 `.tpr`、`.xtc` 以及最重要的 `frame_*.pdb` 文件。
3. **一键运行**：在 Windows 的终端中执行：
   ```bash
   python auto_prolif.py
   ```
4. **收取结果**：终端会实时打印每个体系的处理进度（如 `[Hit1] ✅ 处理成功！`）。全部跑完后，直接去脚本同级目录下的对应靶点文件夹里查看图表和 CSV 数据即可。

***

## 📈 工具三：`auto_rmsf.py` (RMSF 均方根涨落提取与绘图)

### 1. 工具简介
这个脚本用来批量计算复合物体系的 RMSF（均方根涨落），以直观反映配体结合前后，蛋白特定区域（如柔性 Loop 区或催化口袋）的刚性化程度与运动学变化。

**它的核心作用是打通 WSL 终端与 Windows 的数据处理壁垒，并提供了灵活的交互选项：**
* **交互式组别选择**：启动脚本时会在终端弹出一个交互菜单，让你自由选择基于哪个组别进行计算：
  * `[1] Protein`：全体蛋白（包含侧链，波动及噪音较大）。
  * `[3] C-alpha`：碳阿尔法（去除侧链干扰，最能反映骨架真实运动，推荐用于发表文章）。
  * `[4] Backbone`：主链骨架（平滑且稳定，可作为备选）。
* **智能防覆盖命名**：根据你的选择，脚本会自动给输出文件附加相应的后缀（如 `_RMSF_calpha.csv`）。这意味着你可以对同一个体系跑多次不同组别的计算，文件不会互相覆盖。
* **数据提纯与可视化**：自动剥离 GROMACS 原生 `.xvg` 文件中的冗余注释，提取出纯净的二维坐标存为 `.csv` 文件，并利用 Python 绘制带有阴影填充区域的学术级折线图（`.png`），彻底抛弃简陋的默认 Grace 图表。

### 2. 目录结构与依赖文件
程序运行时，会读取 `config.yaml` 里的 `base_dir_win` 目录，并在该目录下自动建立一个专门的 `RMSF` 文件夹存放所有结果。

请确保你的原始数据按照以下结构存放：

```text
E:\GROMACS\Two Target\          <-- 你的工程根目录 (config.yaml 中的 base_dir_win)
├── STAT3\                      <-- 靶点名称 (config 中的 targets)
│   ├── Hit1\                   <-- 具体的复合物体系 (以 config 中的 hit_prefix 开头)
│   │   ├── md_fit.xtc          <-- 【必须】已去除 PBC 并居中叠合的轨迹文件
│   │   └── md.tpr              <-- 【必须】MD 运行输入拓扑文件
│   └── Hit2\
└── RMSF\                       <-- 【运行后自动生成】存放所有RMSF结果的总目录
    └── STAT3\
        ├── Hit1\
        │   ├── rmsf_calpha.xvg           <-- 原始输出的 XVG 文件
        │   ├── Hit1_RMSF_calpha.csv      <-- 提取出的干净数据表
        │   └── Hit1_RMSF_calpha.png      <-- 渲染好的折线图
        └── Hit2\
```

### 3. 操作指南
1. **核对配置**：打开 `config.yaml`，确保 `rmsf` 节点下的 `conda_env`（能调用 `gmx` 命令的环境）以及 `group_default`（如果你不想每次都手动选，可以直接回车使用此默认组别）设置正确。
2. **检查文件**：确认每个待计算的 Hit 目录中都包含 `.xtc` 轨迹文件和 `.tpr` 文件。
3. **一键运行**：在 Windows 终端中执行：
   ```bash
   python auto_rmsf.py
   ```
4. **交互选择**：终端会提示你输入要计算的组别编号（输入 `1`、`3` 或 `4`，直接敲击回车则使用 yaml 文件中的默认值）。
5. **收取结果**：脚本会自动调用 WSL 跑完所有体系。完成后，去工程根目录下新生成的 `RMSF` 文件夹中，按靶点和 Hit 名称查看对应的图表与数据源。

***

## 🧬 工具四：`auto_hbonds.py` (动态氢键数量提取与绘图)

### 1. 工具简介
分子间氢键的数量是衡量配体与蛋白结合稳定性的重要指标。本工具使用 `gmx hbond-legacy` 命令，全自动提取整条分子动力学轨迹中配体（Ligand）与受体（Protein）之间形成的氢键数量。

**核心特性：**
* **零配置自动读取**：脚本极其聪明地复用了 `config.yaml` 中为 `mmpbsa` 设置的 `group_receptor` 和 `group_ligand` 编号，自动模拟键盘输入选择蛋白和配体，免去额外配置的烦恼。
* **单位自动换算**：在绘制图表时，自动将 GROMACS 默认的皮秒（ps）时间轴换算为学术界发文更常用的纳秒（ns）。
* **高级质感出图**：抛弃了简陋的线条，针对氢键这种高频波动的数据，采用了带有暖色调（红色系）半透明填充区域的学术风折线图，使其与冷色调的 RMSF 图表形成完美的视觉区分。

### 2. 目录结构与依赖文件
程序运行时，会在工程根目录下自动建立一个专门的 `HBonds` 文件夹存放所有结果。

请确保原始数据中包含 `.xtc`、`.tpr` 以及 `.ndx` 索引文件：

```text
E:\GROMACS\Two Target\          <-- config.yaml 中的 base_dir_win               
├── STAT3\                      
│   ├── Hit1\                   
│   │   ├── md_fit.xtc          <-- 【必须】已叠合好的轨迹文件
│   │   ├── md.tpr              <-- 【必须】MD 运行输入文件
│   │   └── index.ndx           <-- 【必须】包含受体和配体组别的索引文件
│   └── Hit2\
└── HBonds\                     <-- 【运行后自动生成】存放氢键结果的总目录
    └── STAT3\
        ├── Hit1\
        │   ├── total_hbonds.xvg          <-- 原始数量文件
        │   ├── lifetime.xvg              <-- 寿命文件 (仅生成备用)
        │   ├── Hit1_HBonds.csv           <-- 提取出的时间-数量干净数据表
        │   └── Hit1_HBonds.png           <-- 渲染好的氢键波动图
        └── Hit2\
```

### 3. 操作指南
1. **核对配置**：只要你之前在 `config.yaml` 填好了 `targets`、`hit_prefix` 以及 `mmpbsa` 下面的受体配体组别号，这步直接跳过。
2. **一键运行**：在 Windows 终端中执行：
   ```bash
   python auto_hbonds.py
   ```
3. **收取结果**：脚本跑完后，去根目录新生成的 `HBonds` 文件夹中，按靶点和 Hit 名称收取 `.xvg` 原始文件、供 Origin 使用的 `.csv` 数据源以及高清的 `.png` 分析图。

***

# 工具五 decomp_get.py 结合能单残基能量分解与最强氨基酸统计

## 1. 工具简介
该脚本用于深入分析 MMPBSA 计算产生的 FINAL_DECOMP_MMPBSA.dat 文件，提取配体与受体相互作用中的单残基能量贡献。脚本会自动剥离复合物和受体的绝对能量，精准定位到 DELTAS 段落，提取平均结合能及其标准误。程序会自动将原始的残基命名转换为标准格式如 GLN843。除了为每个体系生成前6个关键氨基酸的柱状图和对应的数据表外，脚本新增了全局统计功能，自动将所有体系中结合能最强的那个氨基酸提取出来，汇总成一张总表，方便后续统一制表对比。

## 2. 目录结构与依赖文件
程序运行时会读取工程根目录下 MMPBSACal 文件夹中的计算结果，并将清洗后的分析数据统一输出到 Dresult 文件夹中。

```text
E:\GROMACS\Two Target\
├── Tool\
│   └── decomp_get.py
├── MMPBSACal\
│   └── STAT3\
│       └── Hit1\
│           └── FINAL_DECOMP_MMPBSA.dat
└── Dresult\
    ├── Strongest_Residue_Summary.csv
    └── STAT3\
        └── Hit1\
            ├── Hit1_All_Residues_Decomp.csv
            ├── Hit1_Top6_Residues.csv
            └── Hit1_Top6_Energy.png
```

## 3. 操作指南
在 Windows 终端中执行 python Tool/decomp_get.py 即可。程序会自动遍历所有靶点和体系，提取并绘制能量柱状图。所有单独体系的图表和数据均存放在 Dresult 对应的文件夹下，用于跨体系对比的单氨基酸最强结合能汇总表则直接生成在 Dresult 的根目录中。

# 工具六 auto_fel.py PCA-FEL 发表级自由能形貌图一键生成工具

## 1. 工具简介
本脚本专为生成符合高质量 SCI 论文标准的自由能形貌图（FEL）而设计。
与传统依赖配体 RMSD 和 Rg 这种容易产生剧烈波动的局部参数不同，本脚本采用**主成分分析（PCA）**策略：通过提取受体蛋白骨架（Backbone）的前两个主成分（PC1 和 PC2）作为坐标轴，能够完美滤除局部噪音，展现复合物在模拟过程中的核心构象转变与真实能量谷。

**核心特性：**
* **严格的 10 步标准化流程**：全程自动处理周期性边界（PBC）、消除平移旋转、剥离复合物，并锁定 Backbone 组别进行协方差矩阵计算，彻底告别配体“瞬移”导致的图像伪影。
* **“混合动力”执行架构**：巧妙结合双系统优势。前端在 WSL (Linux) 环境下调用 GROMACS 算力处理海量轨迹数据；后端自动将结果传回 Windows 环境，利用成熟的 Conda 绘图环境（Matplotlib）直接渲染出图，完美避开 WSL 缺包报错。
* **实时流式监控**：脚本采用了流式输出设计，在 Windows 终端中可实时查看 GROMACS 处理每一帧的进度百分比，告别盲目等待。
* **断点续传**：支持智能识别已生成的体系，中断后再次运行会自动跳过已完成的 Hit。

## 2. 目录结构与依赖文件
该工具高度依赖 `Script` 目录下的三个辅助绘图脚本。运行前请确保目录结构如下，并且每个体系的 `index.ndx` 索引文件中包含完整的 `[ Protein_Lig ]` 组（通常为第 20 组）。

```text
E:\GROMACS\Two Target\
├── config.yaml                 <-- 全局配置文件
├── Script\                     <-- 【必须】辅助脚本目录
│   ├── pc_combine.py           <-- 合并主成分数据
│   ├── xpm2all.bsh             <-- 转换 3D 数据格式 (Linux 脚本)
│   └── xpm2png.py              <-- 无窗口版 Matplotlib 渲染脚本
├── Tool\
│   └── auto_fel.py             <-- 本执行脚本
├── STAT3\                      
│   └── Hit1\                   
│       ├── md.xtc              <-- 【输入】原始轨迹
│       ├── md.tpr              <-- 【输入】运行参数
│       └── index.ndx           <-- 【输入】索引文件 (需包含 Protein_Lig 组)
└── FEL\                        <-- 【运行后自动生成】
    └── STAT3\
        └── Hit1\
            ├── pc1.xvg         <-- 第 1 主成分投影
            ├── pc2.xvg         <-- 第 2 主成分投影
            ├── pc12_gibbs.xpm  <-- 吉布斯自由能矩阵
            ├── pc12_gibbs.xyz  <-- 用于 Origin 等第三方软件的 3D 数据
            └── pc12_gibbs.png  <-- 最终发表级 FEL 渲染图
```

## 3. 操作指南
在 Windows 终端（如 Anaconda Prompt 或 PowerShell）中，激活你的主环境（如 `prolif_env`）后，直接执行以下命令：

```bash
python Tool/auto_fel.py
```

**运行过程说明：**
1. 脚本读取 `config.yaml` 中的靶点列表，自动向 WSL 发送指令。
2. 屏幕会实时打印 `[1/10]` 到 `[9/10]` 的 GROMACS 运行进度（提取复合物 -> 拟合 -> 算协方差 -> 算自由能）。
3. 当单个体系的计算完成后，脚本会立刻在 Windows 侧接管，打印 `🎨 [Windows 绘图] 正在渲染 PNG 图片...`。
4. 所有生成的中间数据和图片，均会按 `靶点/Hit` 的分类，整齐地输出在根目录下的 `FEL` 文件夹中。