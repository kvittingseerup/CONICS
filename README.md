# phyngle
__*PHY*__logenies and clone-specific expression from si__*NGLE*__-cell RNA sequencing

## Table of contents
- [Identifying CNVs from scRNA-seq](#Calling_CNV)
- [Phylogenetic tree contruction](#Constructing_Tree)
- [Intra-clone co-expression networks](#CX_Net)
- [Assessing the correlation of CNV status with single-cell expression](#Corr)
- [False discovery rate estimation: Cross validation](#10x)
- [False discovery rate estimation: Empirical testing](#Empirical)


## <a id="Calling_CNV"></a> Identifying CNVs from scRNA-seq 

### Requirements
  * [Python](https://www.python.org) and [Perl](https://www.perl.org)
  * [beanplot](https://www.jstatsoft.org/article/view/v028c01)(R package)
  * [samtools](http://www.htslib.org)
  * [bedtools](http://bedtools.readthedocs.io/en/latest)  
  * Two directories, the first containing the aligned scRNA-seq data to be classified by CNV status, and a second, containing aligned scRNA-seq data to be used as a control.
  * A file contianing the genomic coordinates of the CNVs in [BED](https://genome.ucsc.edu/FAQ/FAQformat#format1) format.

### Config file
Adjust __Phyngle.cfg__ to customize the following:
  * Path to python/samtools/bedtools/Rscript
  * Thresholds for mapping-quality and read-count
  * FDR for CNV calling

### Running

  ```
  bash run_Phyngle.sh [directory for tumor] [directory for normal] [.bed file for CNV segments] [base name]
  ```
  * __[directory for tumor]__: path to directory containing aligned bam files to be tested. Example glioblastoma data, used in the manuscript, can be obtained [here](https://www.ebi.ac.uk/ega/studies/EGAS00001002185).
    
  * __[directory for normal]__: path to directory containing aligned bam files to be used as a control. Example nonmalignant brain data, used in the manuscript, can be obtained [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835) was used as an examples for the journal 
   
  * __[.bed file for CNV segments]__: tab-delimited bed file of CNV segments to be quantified.
  
      ```
      [chromosome]	[start]	[end]	[chromosome:start:end:CNV]
      ```
      
      __Note: the 4th column of the file must have the exact format shown here:
        __ (__Amp__:amplification, __Del__:deletion)
    * example (SF10281c.cnv.merged_gt1500000_20percent.bed)
  ```
  7   19533   157408385   7:19533:157408385:Amp
  9   19116859    32405639    9:19116859:32405639:Del
  ```
  * __[base name]__ : base name for output directory
  

### Output
All output files will be located in the directory __output_[base name]__.
  1. __incidenceMatrix.csv__: matrix of presence/absence for all CNVs, in individual cells
  2. __pdf__: 
    1. Read-count distribution in CNV segments. (violin plot)
    2. Hierarchical clustering of the single cells by CNV status.

![violin](images/Phyngle_violin.jpg?raw=true "violin" )
![dendrogram](images/Phyngle_dendrogram.jpg?raw=true "dendrogram" )


## <a id="Constructing_Tree"></a> Phylogenetic tree contruction
Phyngle can generate a phylogenetic tree from the CNV incidence matrix, using the Fitch-Margoliash algorithm. Other phylogenetic reconstruction algorithms can be applied, using the incidence matrix as a starting point.

### Requirements
  * [Rscript](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html)
  * [Rphylip](https://cran.r-project.org/web/packages/Rphylip/index.html) (R package)
  * [Phylip](http://evolution.genetics.washington.edu/phylip.html) 

### Config file

Adjust __Tree.cfg__ to change the following.
  * Path to Rscript
  * Path to Rphylip
  
### Running
  * Before running, set the path to Phylip in __Tree.cfg__ file.
```
bash run_Tree.sh [CNV presence/absence matrix][number of genotypes] [base name for output file]
```
  * __[CNV presence/absence matrix]__: .incidenceMatrix.csv files. 
  * __[number of genotypes]__: the number of genotypes to model
  * __[base name]__ : base name for output directory

### Output
__[base name]_cluster.pdf__ (phylogenetic trees) and __[base name]_cluster.txt__ will be generated in the __output_[base name]__ directory. Each leaf corresponds to a clusters of cells with a common genotype. Cluster assignments for each cell will be in __[base name]_cluster.txt__ . 



![tree](images/Trees_cluster.jpg?raw=true "tree" ) 

```
cluster_1  D12,E10,F9,G3,A12,C8,C9,A3,A5,A6,C3,C2,C1,C7,H12,C4,D8,D9,A9,E4,E7,E3,F1,E1,B5,B7,E9,B3,D7,D1
cluster_2  E8,G7,G9,A7,G2,B6,E2
cluster_3  H3,A2,A4,H8,G11,F2,F3,H1,H7
cluster_4  A10,B2
cluster_5  C5
cluster_6  F8,B1
```

## <a id="CX_Net"></a> Intra-clone co-expression networks
Phyngle can construct the local co-expression network of a given gene, based on correlations across single cells.  
  

### Requirements
  * [scde](http://hms-dbmi.github.io/scde)
  * [PCIT](https://cran.r-project.org/web/packages/PCIT/index.html)
  * [boot](https://cran.r-project.org/web/packages/boot/)
  * [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)
  * [raster](https://cran.r-project.org/web/packages/raster/)
  * [flashClust](https://cran.r-project.org/web/packages/flashClust/index.html)
  
### Config file
Adjust __CorrelationNetwork.cfg__ to configure the following:
  * Path to Rscript
  * ncore: Number of cores (default: 12)
  * cor_threshold: Starting threshold to construct the co-expression network (default: 0.9)
  * min_neighbours: How many direct neighbours of gene of interest should be analyzed (default: 20)
  * minRawReads: How many raw reads should map to a gene for it to be included (default: 100)
  * percentCellsExpressing: Percentage (0.15 =15%) of cells expressing a gene for it to be included (default: 0.15)
  * minGenesExpr: How many genes should be expressed in a cell for it to be included (default: 800)
  * depth: How deep should the gene analysis search. (2=only direct neighbor genes would be considered) (default: 2)


### Running

  ```
  bash run_CorrelationNetwork.sh [input matrix] [centered gene] [base name]
  ```
  
  * __[input matrix]__: tab-delimited file of read counts for each gene (rows), for each cell (columns). 
  * __[centered gene]__: a target gene of which neighbor genes are analyzed.  
  * __[base name]__ : base name for output directory

### Output
All the output files will be located in __output_[base name]__.
  1. __[base name]_[correlstion_threshold]_[gene_name].txt__ : co-expression network
  2. __[base name]_[gene_name]_corMat.rd__: Rdata containing the adjusted correlation matrix
  3. __[base name]_topCorrelations.pdf__: bar graph of top correlations. 
  
![CXnet](images/PTEN_topCorr.jpg?raw=true "CXnet" )



## <a id="Corr"></a> Assessing the correlation of CNV status with single-cell expressionq 

### Requirements
  * [zoo](https://cran.r-project.org/web/packages/zoo/index.html)
  
### Config file
Adjust __CompareExomeSeq_vs_ScRNAseq.cfg__ to set the following:
  * Path to Rscript
  * window size for assessing CNV status
  
### Running

  ```
  bash run_compareExomeSeq_vs_ScRNAseq.sh [matrix for read counts] [base name for output file]
  ```
  * __[matrix for read counts]__: tab-delimited file of the number of mapped reads to each gene in the DNA sequencing and in scRNA-seq
  
      ```
      [gene] [chromosome] [start] [#read in DNA-seq(normal)] [#read in DNA-seq(tumor)] [#read in scRNA-seq(normal)] [#read in scRNA-seq(tumor)]
      ```
    * example
  ```
  
DDX11L1   1   11874   538   199   5   0
WASH7P    1   14362   4263   6541   223   45
  ```
  * __[base name]__ : base name for output directory



### Output
__Compare_[window_size].pdf__ (Box plot) will be generated in the directory __output_[base name]__. 
![compare](images/Compare_200.jpg?raw=true "compare" )


## <a id="10x"></a> False discovery rate estimation: Cross validation
Phyngle can estimate false discovery rate via 10-fold cross-validation, using the user-supplied control scRNA-seq dataset. For example, in the manuscript cross validation was performed using [normal brain controls](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835).

### Requirements
  * [beanplot](https://www.jstatsoft.org/article/view/v028c01)(R package)
  * [samtools](http://www.htslib.org)
  * [bedtools](http://bedtools.readthedocs.io/en/latest)  


### Config file
Adjust __10X_cross_validation.cfg__ to set the following:
  * Paths to python/samtools/bedtools/Rscript
  * Thresholds for mapping-quality and read-count.
  * FDR for CNV calling

### Running

  ```
  bash run_10X_cross_validation.sh [directory for control scRNA-seq] [.bed file containing CNV segments] [base name]
  ```
  * __[directory for test]__: path to directory containing the aligned BAM files of the scRNA-seq control data.
   
  * __[.bed file for CNV segments], [base name]__ : same as described in run_Phyngle.sh
;
### Output
Box plot of 10 FDRs resulting from each pooled sample would be generated (__[base name]_boxplot.pdf__) in  __output_[base name]__ directories.
![10X](images/10X_boxplot.jpg?raw=true "10Xval_Test" )

## <a id="Empirical"></a> False discovery rate estimation: Empirical testing
The FDR could be estimated by empirical testing. The number of false positive CNV calls was calculated from non-malignant [fetal brain dataset](http://dx.doi.org/10.1016/j.cell.2015.09.004) which are independent on the [training set](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835)

### Requirements
  * [beanplot](https://www.jstatsoft.org/article/view/v028c01)(R package)
  * [samtools](http://www.htslib.org)
  * [bedtools](http://bedtools.readthedocs.io/en/latest)
  
### Config file
Adjust __Empirical_validation.cfg__ to change the following:
  * Paths to python/samtools/bedtools/Rscript
  * Thresholds for mapping-quality and read count
  * FDR for CNV calling

### Running

  ```
  bash run_empirical_validation.sh [directory for train] [directory for test] [.bed files for CNV segments] [base name]
  ```

  * __[directory for train]__: path to directory containing aligned bam files of scRNA-seq data used as a control to call CNVs
    
  * __[directory for test]__: path to directory containing aligned bam files of scRNA-seq data known not to have CNVs, used as a gold standard.
   
  * __[.bed file for CNV segments], [base name]__ : same as described in run_Phyngle.sh

### Output
Box plot of FDRs will be generated (__[base name]_boxplot.pdf__) in the  __output_[base name]__ directory.
![empirical](images/Empirical_boxplot.jpg?raw=true "empirical_Test" )
