seqfile = <你的序列比对文件名.phy>
treefile = /home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk
outfile = <geneX_AE_bsA_alt_results.txt> * 修改输出文件名

noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
icode = 0
clock = 0
Small_Diff = 1e-5

model = 2       * 指定前景分支有不同的omega类别
NSsites = 2     * 使用允许正选择的位点模型 (M2a的位点部分，应用于前景枝)
fix_omega = 0   * 估计omega
omega = 2.0     * 正选择类别omega的初始值 (修改为2.0)

fix_kappa = 0
kappa = 2
fix_alpha = 1
alpha = .0
Malpha = 0
ncatG = 4
getSE = 0
RateAncestor = 0
method = 0
fix_blength = 2 