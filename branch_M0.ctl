seqfile = <你的序列比对文件名.phy>
treefile = /home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk * 对于M0模型，树上的#1标记会被忽略，但使用同一个树文件是规范的
outfile = <geneX_AE_branch_M0_results.txt> * 修改输出文件名

noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
icode = 0
clock = 0

model = 0       * 所有分支共享一个omega
NSsites = 0     * 所有位点共享一个omega
fix_omega = 0   * 估计omega
omega = 0.4     * 初始omega值

fix_kappa = 0
kappa = 2
fix_alpha = 1
alpha = .0
Malpha = 0
ncatG = 4
getSE = 0
RateAncestor = 0
method = 0
fix_blength = 0