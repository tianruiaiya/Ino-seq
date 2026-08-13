# Ino-seq v1.0.0 中文端到端使用指南

本指南对应一条完整正式流程：从实验样本和 blank/control 的双端 FASTQ
开始，自动完成模块 01–05、队列 QC、脱靶汇总和全流程验收。

正式使用只需要记住根目录的统一命令：

```bash
./inoseq <command>
```

## 1. 先理解一条主线

```text
Phase A：逐样本处理（模块 01–02）
  FASTQ → reads QC → UMI 提取 → 比对 → UMI consensus → ABE signature reads

Phase B：实验/对照配对检测（模块 03–05）
  signature reads → 背景过滤 → 候选区间 → sgRNA 比对 → 位点分类

Phase C：队列汇总与验收
  样本 QC + 脱靶汇总 → TSV/Excel → 完成状态
```

三者是同一个流程的执行阶段，不是三套可任选的分析方案。代码中的
`step01` 对应 Phase A，`postprocess` 对应 Phase B；这些名称仅为兼容旧版
脚本和目录而保留。

## 2. 创建运行环境

在仓库根目录执行：

```bash
mamba env create -f envs/inoseq.yml
conda activate ino-seq
```

已有环境需要同步时：

```bash
mamba env update -n ino-seq -f envs/inoseq.yml --prune
```

## 3. 初始化本地配置

```bash
./inoseq init
```

该命令只在文件不存在时创建以下三个本地配置，不会覆盖已有文件：

```text
config/inoseq.env
config/samples.tsv
config/pairs.tsv
```

三者关系如下：

| 文件 | 内容 | 核心约束 |
|---|---|---|
| `inoseq.env` | 参考基因组、输出目录、资源和参数 | `REFERENCE_FASTA` 必须可读 |
| `samples.tsv` | 所有实验和 control 的 FASTQ | 每个样本 ID 唯一 |
| `pairs.tsv` | 实验样本、匹配 control 和 sgRNA | 所有 ID 必须已在 `samples.tsv` 定义 |

### 3.1 设置运行路径

至少修改：

```bash
REFERENCE_FASTA=/absolute/path/to/hg38.fa
OUTPUT_DIR=/absolute/path/to/new_output
LOG_DIR=/absolute/path/to/new_logs
```

同一次分析失败后继续运行时，应保留原来的 `OUTPUT_DIR`，流程才能识别并复用
已成功阶段；`LOG_DIR` 也建议保持不变以便连续追踪。只有在开始一个科学上独立
的新分析时，才应设置新的 `OUTPUT_DIR` 和 `LOG_DIR`。

### 3.2 填写样本表

`samples.tsv` 必须同时包含实验样本和 blank/control：

```text
sample_id<TAB>read1<TAB>read2
sampleA<TAB>/absolute/path/sampleA_R1.fastq.gz<TAB>/absolute/path/sampleA_R2.fastq.gz
blankA<TAB>/absolute/path/blankA_R1.fastq.gz<TAB>/absolute/path/blankA_R2.fastq.gz
```

### 3.3 填写配对表

```text
sample_id<TAB>control_id<TAB>sgrna
sampleA<TAB>blankA<TAB>ACGTACGTACGTACGTACGTNGG
```

当前版本支持多个实验样本、不同或共用对照，以及不同或共用 sgRNA；同一
实验 `sample_id` 在 `pairs.tsv` 中只能出现一次。同一样本多 sgRNA 或 pooled
sgRNA 文库不属于 v1.0.0 输入模型。

## 4. 准备参考基因组

```bash
./inoseq prepare-reference
```

程序检查并按需准备 BWA 索引、FASTA 索引和序列字典。若参考目录只读，
应由管理员预先提供完整索引。

## 5. 先验证，再查看任务图

完整预检：

```bash
./inoseq validate
```

预检同时核对：软件环境、参考文件、FASTQ、样本表格式、配对关系、sgRNA
格式和关键阈值。

只打印断点续跑后的实际 Slurm 依赖图，不提交任务：

```bash
./inoseq plan
```

输出中的 `[REUSE]` 表示结果通过“完成标记 + 必需输出 + 指纹”检查，将被跳过；
`[RUN]` 表示需要提交。确认每个实验样本只依赖自身及匹配 control，并能看到
所需 QC 和 finalizer 任务后，再正式提交。

## 6. 一次提交完整流程

```bash
./inoseq submit
```

该命令同时适用于第一次运行和失败后的再次提交。默认 `auto` 模式会逐一判断：

```text
每个样本：模块01已完成？模块02已完成？
每个实验/对照：模块03已完成？模块04已完成？模块05已完成？
队列层面：Phase A QC、Phase B QC、finalize 是否仍然有效？
```

已完成且指纹一致的阶段不会重新提交；失败、缺文件或指纹变化的阶段会从最早
无效模块开始，并自动补做其下游汇总和最终验收。

该命令自动建立：

```text
所有样本 Phase A 作业
  ├─→ 样本队列 QC
  └─→ 每个实验/对照的 Phase B 作业
         └─→ 脱靶队列汇总
                 └─→ Phase C 全流程验收
```

所有下游使用 Slurm `afterok` 依赖。任一上游失败时，相关下游不会把不完整
结果标记为成功。

## 7. 查看任务与完成状态

统一状态命令：

```bash
./inoseq status
```

它会显示版本、输出路径、是否完成、最终 Slurm job ID、任务图文件位置，以及
模块01–05、两类队列 QC 和 finalize 的 `CURRENT/INCOMPLETE/STALE` 状态。

使用测试配置或非默认配置时：

```bash
./inoseq status \
  config/samples.test.tsv \
  config/pairs.test.tsv \
  config/inoseq.test.env
```

查看最终作业：

```bash
final_job=$(cat logs/inoseq_last_full_workflow_job_id.txt)
squeue -j "$final_job"
sacct -j "$final_job" \
  --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

查看全部作业依赖：

```bash
column -t -s $'\t' logs/inoseq_last_full_workflow_jobs.tsv
```

正式完成必须同时满足：

```text
最终作业 State=COMPLETED、ExitCode=0:0
<OUTPUT_DIR>/INOSEQ_WORKFLOW_COMPLETE 存在
<OUTPUT_DIR>/QC/full_workflow_status.tsv 中 status=COMPLETED
```

## 8. 从哪里读取结果

优先从 `OUTPUT_DIR/QC/` 读取队列结果：

| 文件 | 用途 |
|---|---|
| `full_workflow_status.tsv` | 全流程状态、版本、样本数、配对数和核心阈值 |
| `inoseq_qc_summary.tsv` | FASTQ、比对、UMI 和 signature-read QC |
| `inoseq_offtarget_summary.tsv` | 各样本 on-target/dependent/independent 位点汇总 |
| `inoseq_strand_summary.tsv` | 链方向和 protospacer 分类汇总 |
| `dependent_target_analysis.xlsx` | 可直接查看的队列 Excel 报告 |

单样本最终位点表：

```text
<OUTPUT_DIR>/<sample>/postprocess/offtarget/<sample>_dependent_target.txt
```

样本实际运行参数：

```text
<OUTPUT_DIR>/<sample>/postprocess/run_parameters.tsv
```

全部中间文件和字段含义见 `docs/OUTPUT_CONTRACT.md`。

## 9. 失败时按任务图定位

1. 在 `logs/inoseq_last_full_workflow_jobs.tsv` 找到失败阶段和 job ID。
2. 用 `sacct` 确认 `State`、`ExitCode` 和资源使用。
3. 查看对应的 `logs/<job_name>_<job_id>.err` 与 `.out`。
4. 检查该阶段完成标记是否缺失，不要仅凭某个中间文件存在就判断成功。
5. 修复错误后再次执行 `./inoseq plan`，确认恢复起点，再执行
   `./inoseq submit`；不需要手动删除全部输出。

| 失败位置 | 首要检查 |
|---|---|
| Phase A | FASTQ、参考索引、fastp/cutadapt/BWA/fgbio、内存 |
| Phase B | 实验与 control 的 `.end`/`_end.bam`、配对表和 sgRNA |
| 队列 QC | 单样本完成标记和必需 QC 文件 |
| finalizer | 所有启用分支的完成标记及队列 TSV/Excel |

## 10. 自动续跑、指定起点与旧结果接纳

### 10.1 默认自动续跑

```bash
./inoseq plan
./inoseq submit
```

这是日常恢复失败任务的首选方式。不同样本可以从不同模块恢复，例如 control
全部复用、实验样本从模块02继续、配对分析从模块03重新开始。

### 10.2 明确指定从哪一步开始

先用 `plan` 审核，再以相同选项提交：

```bash
./inoseq plan   --from-stage module04
./inoseq submit --from-stage module04
```

| 起点 | 行为 |
|---|---|
| `module01` | 所有样本重做模块01–02，随后重做模块03–05、QC和验收 |
| `module02` | 复用有效模块01，重做模块02及全部下游 |
| `module03` | 复用有效Phase A，重做模块03–05及下游 |
| `module04` | 复用有效模块03，重做模块04–05及下游 |
| `module05` | 复用有效模块04，只重做模块05、脱靶汇总和验收 |
| `qc` | 不重做模块01–05，只重做两类队列汇总和验收 |
| `finalize` | 只重新执行最终完整性验收 |

`--force` 等价于 `--from-stage module01`。指定较晚起点时，如果其前置模块
缺失或指纹已变化，程序会拒绝提交并提示应从哪个更早模块开始，不会带病复用。

### 10.3 指纹检查范围

阶段只有同时满足以下三项才显示为 `CURRENT`：

1. 该阶段完成标记存在；
2. 输出契约规定的必需文件全部存在；
3. 当前指纹与成功时记录一致。

指纹包括相关 FASTQ/参考文件元数据、配置参数、上游指纹、模块代码和 Ino-seq
版本。大 FASTQ 与参考基因组使用路径、大小和修改时间，不在每次提交时重新计算
整文件 SHA256；流程代码和小型样本/配对表使用内容 SHA256。

### 10.4 接纳升级前已有的完整结果

旧版 `v1.0.0` 结果只有完成标记、没有模块指纹，因此不会被自动静默复用。确认
这些结果确实对应当前 FASTQ、参考基因组和配置后，只执行一次：

```bash
./inoseq adopt-existing
./inoseq status
./inoseq plan
```

`adopt-existing` 不运行任何分析模块，只检查必需输出并生成状态记录。如果输出
不完整，它会失败且不会把缺失阶段伪装为成功。

## 11. 冻结参数与测试

v1.0.0 核心筛选条件：

```text
切割位点窗口：±15 bp
fold change：≥1.5
原始 P 值：<0.05
候选位点合并距离：≤30 bp
候选最小总切割 reads：≥3
sgRNA 搜索窗口：±25 bp
最大 sgRNA 比对得分：8
邻近位点扫描：±100 bp
```

BH-FDR 仅输出、不参与过滤；没有额外的“实验组 reads >5”过滤条件。

运行仓库测试：

```bash
./inoseq test
python -m ruff check workflow tests
```

科学计算定义见 `docs/ALGORITHM_CONTRACT.md`，所有输出定义见
`docs/OUTPUT_CONTRACT.md`。
