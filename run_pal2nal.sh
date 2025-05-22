#!/bin/bash

# 脚本功能：批量使用pal2nal.pl将比对后的蛋白质序列和对应的CDS序列转换为密码子比对
# 用法：./run_pal2nal.sh <蛋白质比对目录> <CDS序列目录> <输出目录> [密码子表编号]

# 检查参数
if [ $# -lt 3 ]; then
    echo "用法: $0 <蛋白质比对目录> <CDS序列目录> <输出目录> [密码子表编号]"
    echo "密码子表编号默认为1(通用密码表)，可选值参见pal2nal.pl帮助信息"
    exit 1
fi

# 获取输入参数
ALIGNED_AA_DIR="$1"
CDS_DIR="$2"
OUTPUT_DIR="$3"
CODON_TABLE="${4:-1}"  # 默认使用通用密码表(1)

# 检查目录是否存在
if [ ! -d "$ALIGNED_AA_DIR" ]; then
    echo "错误: 蛋白质比对目录 '$ALIGNED_AA_DIR' 不存在"
    exit 1
fi

if [ ! -d "$CDS_DIR" ]; then
    echo "错误: CDS序列目录 '$CDS_DIR' 不存在"
    exit 1
fi

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 设置pal2nal.pl脚本路径
PAL2NAL="./pal2nal.v14/pal2nal.pl"

# 检查pal2nal.pl是否存在
if [ ! -f "$PAL2NAL" ]; then
    echo "错误: pal2nal.pl脚本不存在于 '$PAL2NAL'"
    echo "请确认pal2nal.pl脚本位置或修改脚本中的PAL2NAL变量"
    exit 1
fi

# 确保pal2nal.pl可执行
chmod +x "$PAL2NAL"

# 计数器
total_files=0
processed_files=0
failed_files=0

# 获取总文件数
total_files=$(find "$ALIGNED_AA_DIR" -type f | wc -l)
echo "开始处理 $total_files 个文件..."

# 循环处理每个蛋白质比对文件
for aa_file in "$ALIGNED_AA_DIR"/*; do
    if [ -f "$aa_file" ]; then
        # 获取文件名（不含路径和扩展名）
        filename=$(basename "$aa_file")
        base_name="${filename%.*}"
        
        # 构建对应的CDS文件名和输出文件名
        cds_file="$CDS_DIR/${base_name}_cds.fa"
        output_file="$OUTPUT_DIR/${base_name}_codon.aln"
        
        echo "处理文件: $filename"
        
        # 检查CDS文件是否存在
        if [ ! -f "$cds_file" ]; then
            echo "  警告: 找不到对应的CDS文件 '$cds_file'，跳过此文件"
            ((failed_files++))
            continue
        fi
        
        # 运行pal2nal.pl
        perl "$PAL2NAL" "$aa_file" "$cds_file" -output paml -codontable "$CODON_TABLE" > "$output_file"
        
        # 检查是否成功
        if [ $? -eq 0 ]; then
            echo "  转换成功: $output_file"
            ((processed_files++))
        else
            echo "  错误: 转换失败"
            ((failed_files++))
        fi
    fi
done

# 输出处理结果
echo "处理完成!"
echo "总文件数: $total_files"
echo "成功处理: $processed_files"
echo "处理失败: $failed_files" 