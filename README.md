# 项目脚本说明

本项目包含一系列用于直系同源基因分析的脚本，从候选基因列表开始，到最终生成比对后的蛋白质序列和转换的codon序列。

## 目录
1.  [run_OrthoFInder_for_66_species.sh]
1.  [1.Search_for_OG_for_candidate_genes.py](#1search_for_og_for_candidate_genespy)
2.  [2. 4_filter_OGs_by_coverage_extract_CDS.py](#2-4_filter_ogs_by_coverage_extract_cdspy)
3.  [3.Translate_cds_to_protein.py](#3translate_cds_to_proteinpy)
4.  [4.run_MAFFT_for_aa_msa.sh](#4run_mafft_for_aa_msash)
5.  [5.pal2nal.pl](#5pal2nalpl)
6.  [6.run_pal2nal.sh](#6run_pal2nalsh)

---

## 1.Search_for_OG_for_candidate_genes.py

**功能**:
此脚本根据提供的候选基因列表和指定的物种，在 OrthoFinder 的输出结果中查找对应的直系同源组 (Orthogroups)，并将这些组相关的多序列比对 (MSA) 文件复制到指定输出目录。

**参数**:
*   `--gene_list_file`: 候选基因列表的 CSV 文件路径。
    *   默认值: `DEG.csv`
*   `--orthogroups_file`: OrthoFinder 输出的 `Orthogroups.tsv` 文件路径。
    *   默认值: `/home/yuhangjia/data/AlternativeSplicing/exon_expansion_test/orthofinder/output_msa_2/Results_Jul05/Orthogroups/Orthogroups.tsv`
*   `--msa_dir`: 包含 MSA 文件的目录路径 (通常是 OrthoFinder 结果中的 `MultipleSequenceAlignments` 目录)。
    *   默认值: `/home/yuhangjia/data/AlternativeSplicing/exon_expansion_test/orthofinder/output_msa_2/Results_Jul05/MultipleSequenceAlignments`
*   `--output_dir`: 保存筛选出的 MSA 文件的目录路径。
    *   默认值: `selected_msa_files_deg`
*   `--species`: 在 `Orthogroups.tsv` 文件中要搜索基因的物种列的名称。
    *   默认值: `Apis_mellifera`

**输入**:
*   一个 CSV 文件，每行包含一个基因 ID（位于第一列）。
*   OrthoFinder 生成的 `Orthogroups.tsv` 文件，其表头包含各个物种的名称。
*   一个包含 OrthoFinder 生成的 MSA 文件的目录。

**输出**:
*   在指定的输出目录 (`--output_dir`) 下，生成一个包含与候选基因相关的直系同源组的 MSA 文件副本。

---

## 3.2._filter_OGs_by_coverage_extract_CDS.py

**功能**:
此脚本用于筛选直系同源组 (Orthogroups) 并提取编码序列 (CDS)。它首先根据用户定义的单拷贝物种覆盖率阈值来筛选在脚本1输出中的 Orthogroups。具体步骤如下：
1.  从 `Orthogroups.tsv` 文件中读取每个 Orthogroup 在各个物种中的基因信息。
2.  排除用户通过 `--excluded_species` 指定的物种 (例如，果蝇)。
3.  基于剩余的物种，计算每个 Orthogroup 中拥有单拷贝基因的物种所占的百分比。
4.  如果此百分比超过 `--single_copy_threshold_percentage` (默认为 70%)，则该 Orthogroup 被视为合格。
5.  对于合格的 Orthogroup，脚本仅保留那些具有单拷贝基因的物种中的基因 (即，如果一个物种在该 Orthogroup 中有多个基因，则这些基因不会被用于提取 CDS)。
6.  最后，为这些筛选出的单拷贝基因提取最长的 CDS 序列。
此脚本有效地取代了旧流程中 `2.Filter_SingleCopyOrthoGroups.py` 和 `3.Extract_longest_cds_for_SCOG.py` 的组合功能，提供了一种基于覆盖率的更精细筛选方法。

**参数**:
*   `--orthogroups_file`: OrthoFinder 输出的 `Orthogroups.tsv` 文件路径。(必需)
*   `--msa_dir_from_script1`: 包含由 `1.Search_for_OG_for_candidate_genes.py` 筛选出的 MSA 文件的目录路径 (例如, `selected_msa_files_deg`)。(必需)
*   `--nucleotide_db_dir`: 包含各物种的核酸序列 (CDS) FASTA 文件的目录路径。每个物种一个文件，文件名即物种名 (例如 `Apis_mellifera.cds.fa`)，FASTA 文件中的序列头部应为 `>gene_id`。(必需)
*   `--output_dir`: 保存提取出的 CDS 序列的目录路径。
    *   默认值: `cds_filtered_by_coverage`
*   `--excluded_species`: 在计算覆盖率和提取 CDS 时需要排除的物种列表 (名称需与 `Orthogroups.tsv` 文件头中的一致)。
    *   默认值: `['Drosophila_melanogaster']` (以空格分隔的列表)
*   `--single_copy_threshold_percentage`: 单拷贝物种覆盖率的最小百分比阈值。
    *   默认值: `70.0`
*   `--subspecies_included`: 一个或多个包含亚种名称的物种全名列表，用于正确解析 MSA 文件头中的基因 ID。
    *   默认值: `['Bombus_vancouverensis_nearcticus', 'Osmia_bicornis_bicornis']` (以空格分隔的列表)

**输入**:
*   OrthoFinder 生成的 `Orthogroups.tsv` 文件。
*   一个包含 MSA 文件的目录 (脚本1的输出，例如 `selected_msa_files_deg`)。
*   一个包含每个物种的 CDS 序列的目录 (FASTA 格式)。

**输出**:
*   在指定的输出目录 (`--output_dir`) 下，为每个通过筛选的 Orthogroup 生成一个对应的 CDS 序列 FASTA 文件。文件名格式为 `[orthogroup_id]_cds.fa`。这些文件仅包含在合格 Orthogroup 中具有单拷贝基因的物种的最长 CDS。

---

## 3.Translate_cds_to_protein.py

**功能**:
此脚本将提取出的 CDS 序列翻译成蛋白质序列。

**参数**:
*   `--input_folder`: 包含 CDS FASTA 文件的目录路径 (脚本2的输出)。
    *   默认值: `cds_filtered_by_coverage`
*   `--output_folder`: 保存翻译后的蛋白质 FASTA 文件的目录路径。
    *   默认值: `translated_proteins_deg`
*   `--genetic_code`: 用于翻译的遗传密码表编号 (NCBI 定义)。
    *   默认值: `1` (标准遗传密码)

**输入**:
*   一个包含 CDS 序列的目录 (FASTA 格式)。

**输出**:
*   在指定的输出目录 (`--output_folder`) 下，为每个输入的 CDS 文件生成一个对应的蛋白质序列 FASTA 文件。文件名中的 `.fa` 或 `.fasta` 会被替换为 `_protein.fa` 或 `_protein.fasta`。

---

## 4.run_MAFFT_for_aa_msa.sh

**功能**:
这是一个 Shell 脚本，使用 MAFFT 工具对翻译后的蛋白质序列进行多序列比对。

**参数 (硬编码在脚本中)**:
*   输入目录: `./translated_proteins_deg` (脚本3的输出)
*   输出目录: `./aligned_translated_proteins_deg`
*   MAFFT 使用的线程数: `64`
*   MAFFT 执行路径: `/usr/local/biotools/m/mafft:7.525--h031d066_0` (通过 Singularity 执行)

**输入**:
*   一个目录 (`./translated_proteins_deg`)，其中包含蛋白质序列的 FASTA 文件 (每个文件通常代表一个直系同源组)。

**输出**:
*   一个目录 (`./aligned_translated_proteins_deg`)，其中包含比对后的蛋白质序列的 FASTA 文件。输出文件名格式为 `[input_filename_prefix]_aligned.fa`。

---

## 5.pal2nal.pl

**功能**:
`pal2nal.pl` (版本 v14) 是一个 Perl 脚本，用于将蛋白质序列比对转换为相应的密码子比对。它能够处理 CLUSTAL 和 FASTA 格式的蛋白质比对，并接受一个或多个包含相应核酸序列的 FASTA 文件。

**用法**:
```bash
perl pal2nal.pl pep.aln nuc.fasta [nuc.fasta...] [options] > output_codon_alignment.aln
```

**参数**:
*   `pep.aln`: 蛋白质比对文件 (CLUSTAL 或 FASTA 格式)。
*   `nuc.fasta`: 一个或多个包含核酸序列的 FASTA 文件 (可以是单个多序列 FASTA 文件或多个单独的文件)。
*   `options`:
    *   `-h`: 显示帮助信息。
    *   `-output <format>`: 指定输出格式。可选格式包括：
        *   `clustal` (默认)
        *   `paml`
        *   `fasta`
        *   `codon`
    *   `-blockonly`: 仅显示在 CLUSTAL 比对中使用 `#` 标记的用户指定区块。
    *   `-nogap`: 移除包含空位和框内终止密码子的列。
    *   `-nomismatch`: 从输出中移除蛋白质和对应 cDNA 之间不匹配的密码子。
    *   `-codontable <N>`: 指定密码子表编号。默认为 `1` (标准通用密码子表)。其他常用选项包括 `2` (脊椎动物线粒体), `5` (无脊椎动物线粒体), `11` (细菌、古菌和植物质体) 等。详细列表请运行 `perl pal2nal.pl -h` 查看。
    *   `-html`: HTML 格式输出 (主要用于 Web 服务器)。
    *   `-nostderr`: 不向 STDERR 输出消息 (主要用于 Web 服务器)。

**输入**:
*   一个蛋白质序列比对文件 (CLUSTAL 或 FASTA 格式)。
*   一个或多个包含与蛋白质序列对应的核酸序列的 FASTA 文件。
    *   **注意**: 输入的蛋白质序列和核酸序列的 ID 应当对应，或者它们的顺序必须在 `pep.aln` 和 `nuc.fasta` 文件中保持一致。脚本会优先使用蛋白质比对文件中的 ID。

**输出**:
*   转换后的密码子比对结果将输出到标准输出 (STDOUT)，可以重定向到文件。输出格式由 `-output` 参数决定。

---

## 6.run_pal2nal.sh

**功能**:
这是一个 Shell 脚本，用于批量使用 pal2nal.pl 将比对后的蛋白质序列和对应的 CDS 序列转换为密码子比对。脚本会自动查找蛋白质比对目录中的所有文件，并在 CDS 目录中寻找对应的核酸序列文件，然后调用 pal2nal.pl 进行转换。

**用法**:
```bash
./run_pal2nal.sh <蛋白质比对目录> <CDS序列目录> <输出目录> [密码子表编号] [-format <输出格式>]
```

**参数**:
*   `<蛋白质比对目录>`: 包含比对后的蛋白质序列的目录路径（如 `aligned_translated_proteins_deg`，脚本4的输出）。
*   `<CDS序列目录>`: 包含 CDS 序列的目录路径（如 `cds_filtered_by_coverage`，脚本2的输出）。
*   `<输出目录>`: 保存生成的密码子比对文件的目录路径。
*   `[密码子表编号]`: 可选参数，指定使用的密码子表编号（默认为 `1`，即标准通用密码子表）。
*   `[-format <输出格式>]`: 可选参数，指定输出格式（默认为 `paml`）。可选值包括：`clustal`、`paml`、`fasta`、`codon`。

**输入**:
*   一个包含比对后的蛋白质序列的目录，其中每个文件通常代表一个直系同源组。
*   一个包含对应 CDS 序列的目录，其中文件名格式应为 `[base_name]_cds.fa`，其中 `[base_name]` 是蛋白质比对文件去掉扩展名后的名称。

**输出**:
*   在指定的输出目录下，为每个成功处理的输入文件生成一个对应的密码子比对文件。输出文件名格式为 `[base_name]_codon.aln`。
*   脚本会显示处理进度和统计信息，包括总文件数、成功处理的文件数和失败的文件数。

**注意事项**:
*   脚本假设 pal2nal.pl 位于 `./pal2nal.v14/pal2nal.pl`。如果您的路径不同，请编辑脚本中的 `PAL2NAL` 变量。
*   请确保您的系统已安装 Perl 环境。

---

## 引文

Mikita Suyama, David Torrents, and Peer Bork (2006)
PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments.
Nucleic Acids Res. 34, W609-W612. 