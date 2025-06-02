seqfile = <你的序列比对文件名.phy>  * 例如: gene1_alignment.phy
treefile = /home/kosukesano/tools/for_paml/IQTREE_6sp/data/new_tree_IQTREE_ultrametric.nwk * 这是你已标记前景枝的树文件
outfile = <具体的输出文件名.txt>    * 例如: gene1_bsA_alt_results.txt

noisy = 9       * 屏幕输出非常详细的信息
verbose = 1     * 输出文件包含详细信息
runmode = 0     * 使用用户提供的树进行ML分析

seqtype = 1     * 密码子序列
CodonFreq = 2   * 使用F3X4模型估算密码子频率 (GY94)
                * (注: 你提供的PAML指南中推荐CodonFreq = 7 (FmutSel) [cite: 122]，但此处遵循你的设置)
icode = 0       * 通用遗传密码

clock = 0       * 不假设分子钟 (各分支进化速率可以不同)
fix_kappa = 0   * 估计kappa (转换/颠换速率比)
kappa = 2       * kappa的初始值
fix_alpha = 1   * 固定alpha值 (不使用gamma分布来模拟位点间速率变异)
alpha = 0.0     * alpha固定为0 (因为fix_alpha=1)
Malpha = 0      * (当alpha=0时，此参数通常不活跃)
ncatG = 4       * (当alpha=0时，此参数通常不活跃)

getSE = 0           * 不计算参数的标准误
RateAncestor = 0    * 不进行祖先序列重建
method = 0          * ML方法中的优化算法 (0: 同时优化所有分支)
fix_blength = 0     * 让codeml在当前模型下估计分支长度
