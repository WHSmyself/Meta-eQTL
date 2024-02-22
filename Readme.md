Meta-eQTL is a new method for eQTL detection. It clusters genes within tissues based on gene expression levels, grouping genes with high correlation together, and performs a Meta-analysis. The method's reliability and practicality were confirmed through analysis of the GTEx.v8 dataset.
1.Code documentation
    The clustering script is written in R language and primarily accomplishes three functions: 1) It clusters genes within the same tissue based on the correlation coefficient of their expression levels; 2) It selects genes of interest and extracts other genes that are clustered with the selected gene in the same category.
    Rscript gene_cluster.R -i input.file -g target_gene -t 0.5 -e 1 -o out_path
    Rscript gene_cluster.R -h
    