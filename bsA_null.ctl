seqfile = <你的序列比对文件名.phy>
treefile = /home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk
outfile = <geneX_AE_bsA_null_results.txt> * 修改输出文件名

noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
icode = 0
clock = 0

model = 2
NSsites = 2
fix_omega = 1   * 固定omega
omega = 1       * 将前景分支上受选择位点的omega固定为1

fix_kappa = 0
kappa = 2
fix_alpha = 1
alpha = .0
Malpha = 0
ncatG = 8
getSE = 0
RateAncestor = 0
method = 0
fix_blength = 2