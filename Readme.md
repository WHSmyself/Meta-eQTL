### Meta-eQTL is a new method for eQTL detection. It clusters genes within tissues based on gene expression levels, grouping genes with high correlation together, and performs a Meta-analysis. The method's reliability and practicality were confirmed through analysis of the GTEx.v8 dataset.

## 1. Code documentation
    
    The clustering script is written in R language and primarily accomplishes three functions: 
    1) It clusters genes within the same tissue based on the correlation coefficient of their expression levels; 
    2) It selects genes of interest and extracts other genes that are clustered with the selected gene in the same category.
    
    Rscript gene_cluster.R -i input.file -g target_gene -t 0.5 -e 1 -o out_path
   ![image](https://github.com/WHSmyself/Meta-eQTL/assets/43985955/c668b2ab-5af2-4c01-a83a-c175b66c3425)
## 2. eQTL analysis

   The eQTL analysis is performed using the GEMMA software, and we have already provided the installation package, which can be obtained from the tools folder.

   gemma -g genotype.file -p phenotype.file -a SNP_Annotation.file -k knishp.file -c cov.file -lmm 1 -n 1 -o outpath
## 3. Meta analysis

   Meta-analysis is conducted using the MTAG and GEMMA software. We have already provided the installation packages for both pieces of software, which can be obtained from the tools folder. 
    
