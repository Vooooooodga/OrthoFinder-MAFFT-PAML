#!/bin/bash

# 脚本功能：批量使用pal2nal.pl将比对后的蛋白质序列和对应的CDS序列转换为密码子比对
# 用法：./run_pal2nal.sh <蛋白质比对目录> <CDS序列目录> <输出目录> [密码子表编号] [-format <输出格式>]

# 默认值
CODON_TABLE="1"  # 默认使用通用密码表(1)
OUTPUT_FORMAT="paml"  # 默认输出格式为paml

# 解析命令行参数
if [ $# -lt 3 ]; then
    echo "用法: $0 <蛋白质比对目录> <CDS序列目录> <输出目录> [密码子表编号] [-format <输出格式>]"
    echo "密码子表编号默认为1(通用密码表)，可选值参见pal2nal.pl帮助信息"
    echo "输出格式默认为paml，可选值: clustal, paml, fasta, codon"
    exit 1
fi

# 获取前三个必要参数
ALIGNED_AA_DIR="$1"
CDS_DIR="$2"
OUTPUT_DIR="$3"
shift 3

# 解析剩余可选参数
while [ $# -gt 0 ]; do
    case "$1" in
        -format)
            if [ -z "$2" ] || [[ "$2" == -* ]]; then
                echo "错误: -format 选项需要一个参数"
                exit 1
            fi
            OUTPUT_FORMAT="$2"
            # 验证输出格式
            if [[ "$OUTPUT_FORMAT" != "clustal" && "$OUTPUT_FORMAT" != "paml" && "$OUTPUT_FORMAT" != "fasta" && "$OUTPUT_FORMAT" != "codon" ]]; then
                echo "错误: 无效的输出格式 '$OUTPUT_FORMAT'。可选值: clustal, paml, fasta, codon"
                exit 1
            fi
            shift 2
            ;;
        -*)
            echo "警告: 未知选项 $1，将被忽略"
            shift
            ;;
        *)
            # 假设是密码子表编号
            if [[ "$1" =~ ^[0-9]+$ ]]; then
                CODON_TABLE="$1"
                shift
            else
                echo "错误: 无效的参数 '$1'"
                exit 1
            fi
            ;;
    esac
done

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

# 输出运行参数摘要
echo "运行参数摘要:"
echo "蛋白质比对目录: $ALIGNED_AA_DIR"
echo "CDS序列目录: $CDS_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "密码子表编号: $CODON_TABLE"
echo "输出格式: $OUTPUT_FORMAT"
echo ""

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
        
        # 运行pal2nal.pl，使用用户指定的输出格式
        perl "$PAL2NAL" "$aa_file" "$cds_file" -output "$OUTPUT_FORMAT" -codontable "$CODON_TABLE" > "$output_file"
        
        # 检查是否成功
        if [ $? -eq 0 ]; then
            echo "  转换成功: $output_file (格式: $OUTPUT_FORMAT)"
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