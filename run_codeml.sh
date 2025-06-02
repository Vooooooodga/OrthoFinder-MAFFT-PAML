#$ -S /bin/bash
#$ -cwd
#$ -l gpu

# ディレクトリの設定
input_dir="/home/kosukesano/tools/for_paml/data/CDS_SCO" #SCOのファイルが格納されているディレクトリを指定する。
bsA_dir="/home/kosukesano/tools/for_paml/241009_IQTREE_6sp/bsA" #対立仮説用のディレクトリを指定する
result_dir="$bsA_dir/result" # 上記のディレクトリで結果を出力する/resultディレクトリを予め作成しておき、指定する。
template_ctl="$bsA_dir/template.ctl" # 対立仮説用のコントロールファイルを指定する。

# 出力ディレクトリが存在しない場合は作成
mkdir -p "$result_dir"

# テンプレートの制御ファイルを読み込む
ctl_template=$(cat "$template_ctl")

# ディレクトリ内の_maffted_fixed.fastaファイルを処理
for file in "$input_dir"/*_maffted.fna; do
  if [[ -f "$file" ]]; then
    base_name=$(basename "$file" .fna)
    outfile_path="$result_dir/${base_name}_branch_alt"

    # 一時的な制御ファイルの内容を生成
    ctl_content="${ctl_template//<SEQFILE>/$file}"
    ctl_content="${ctl_content//<OUTFILE>/$outfile_path}"

    # 一時的な制御ファイルを作成
    ctl_path="$bsA_dir/bsA.ctl"
    echo "$ctl_content" > "$ctl_path"

    # PAMLを実行
    singularity exec -e /usr/local/biotools/p/paml:4.9--h779adbc_6 codeml "$ctl_path"

    echo "Processed file: $file, output: $outfile_path"
  fi
done