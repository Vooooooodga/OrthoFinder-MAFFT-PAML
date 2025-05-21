# 项目脚本说明

本项目包含一系列用于直系同源基因分析的脚本，从候选基因列表开始，到最终生成比对后的蛋白质序列和转换的codon序列。

## 目录
1.  [1.Search_for_OG_for_candidate_genes.py](#1search_for_og_for_candidate_genespy)
2.  [2.Filter_SingleCopyOrthoGroups.py](#2filter_singlecopyorthogroupspy)
3.  [3.Extract_longest_cds_for_SCOG.py](#3extract_longest_cds_for_scogpy)
4.  [4.Translate_cds_to_protein.py](#4translate_cds_to_proteinpy)
5.  [5.run_MAFFT_for_aa_msa.sh](#5run_mafft_for_aa_msash)
6.  [pal2nal.pl](#pal2nalpl)

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

## 2.Filter_SingleCopyOrthoGroups.py

**功能**:
此脚本筛选上一步选出的 MSA 文件，仅保留那些对应于单拷贝直系同源组 (Single-Copy Orthogroups) 的文件。

**参数**:
*   `--single_copy_orthogroups_file`: 列出单拷贝直系同源组 ID 的文本文件路径 (通常是 OrthoFinder 结果中的 `Orthogroups_SingleCopyOrthologues.txt`)。
    *   默认值: `/home/yuhangjia/data/AlternativeSplicing/exon_expansion_test/orthofinder/output_msa_2/Results_Jul05/Orthogroups/Orthogroups_SingleCopyOrthologues.txt`
*   `--selected_msa_files_dir`: 包含由 `1.Search_for_OG_for_candidate_genes.py` 筛选出的 MSA 文件的目录路径。
    *   默认值: `./selected_msa_files_deg`
*   `--output_dir`: 保存筛选出的单拷贝 MSA 文件的目录路径。
    *   默认值: `single_copy_selected_msa_files_deg`

**输入**:
*   一个文本文件，每行包含一个单拷贝直系同源组的 ID。
*   一个包含 MSA 文件的目录 (脚本1的输出)。

**输出**:
*   在指定的输出目录 (`--output_dir`) 下，生成一个仅包含单拷贝直系同源组的 MSA 文件副本。

---

## 3.Extract_longest_cds_for_SCOG.py

**功能**:
此脚本为筛选出的单拷贝直系同源组中的每个基因，从提供的核酸序列数据库中提取最长的编码序列 (CDS)。

**参数**:
*   `--protein_msa_dir`: 包含单拷贝蛋白质 MSA 文件的目录路径 (脚本2的输出)。
    *   默认值: `single_copy_selected_msa_files_deg`
*   `--nucleotide_db_dir`: 包含核酸序列 (CDS) FASTA 文件的目录路径。此目录应每个物种一个 FASTA 文件，文件名即物种名 (例如 `SpeciesA.fasta`)，FASTA 文件中的序列头部应包含基因标识符。
    *   默认值: `../longest_cds`
*   `--output_dir`: 保存提取出的 CDS 序列的目录路径。
    *   默认值: `cds_sequences_deg`
*   `--subspecies_included`: 一个或多个包含亚种名称的物种全名列表，用于正确解析基因 ID。
    *   默认值: `['Bombus_vancouverensis_nearcticus', 'Osmia_bicornis_bicornis']` (以空格分隔的列表)

**输入**:
*   一个包含单拷贝蛋白质 MSA 文件的目录 (FASTA 格式)。
*   一个包含每个物种的 CDS 序列的目录 (FASTA 格式)。

**输出**:
*   在指定的输出目录 (`--output_dir`) 下，为每个输入的 MSA 文件生成一个对应的 CDS 序列 FASTA 文件。文件名格式为 `[orthogroup_id]_cds.fa`。

---

## 4.Translate_cds_to_protein.py

**功能**:
此脚本将提取出的 CDS 序列翻译成蛋白质序列。

**参数**:
*   `--input_folder`: 包含 CDS FASTA 文件的目录路径 (脚本3的输出)。
    *   默认值: `cds_sequences_deg`
*   `--output_folder`: 保存翻译后的蛋白质 FASTA 文件的目录路径。
    *   默认值: `translated_proteins_deg`
*   `--genetic_code`: 用于翻译的遗传密码表编号 (NCBI 定义)。
    *   默认值: `1` (标准遗传密码)

**输入**:
*   一个包含 CDS 序列的目录 (FASTA 格式)。

**输出**:
*   在指定的输出目录 (`--output_folder`) 下，为每个输入的 CDS 文件生成一个对应的蛋白质序列 FASTA 文件。文件名中的 `.fa` 或 `.fasta` 会被替换为 `_protein.fa` 或 `_protein.fasta`。

---

## 5.run_MAFFT_for_aa_msa.sh

**功能**:
这是一个 Shell 脚本，使用 MAFFT 工具对翻译后的蛋白质序列进行多序列比对。

**参数 (硬编码在脚本中)**:
*   输入目录: `./translated_proteins_deg` (脚本4的输出)
*   输出目录: `./aligned_translated_proteins_deg`
*   MAFFT 使用的线程数: `64`
*   MAFFT 执行路径: `/usr/local/biotools/m/mafft:7.525--h031d066_0` (通过 Singularity 执行)

**输入**:
*   一个目录 (`./translated_proteins_deg`)，其中包含蛋白质序列的 FASTA 文件 (每个文件通常代表一个直系同源组)。

**输出**:
*   一个目录 (`./aligned_translated_proteins_deg`)，其中包含比对后的蛋白质序列的 FASTA 文件。输出文件名格式为 `[input_filename_prefix]_aligned.fa`。

---

## pal2nal.pl

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

## 引文

Reference:
Mikita Suyama, David Torrents, and Peer Bork (2006)
PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments.
Nucleic Acids Res. 34, W609-W612. 