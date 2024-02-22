#install.packages("argparse")
#install.packages("getopt")
library("getopt")
library("argparse")
spec <- matrix( c("input", "i", 2, "character", "Input value",
"gene", "g", 2, "character", "Target gene",
"td","t",2,"numeric","Cluster threshold value",
"expression", "e", 2, "numeric", "Whether you need to export gene expression",
"output", "o", 2, "character", "Output path",
"help", "h", 0, "logical", "This is Help!"),
byrow=TRUE, ncol=5) 
opt <- getopt(spec=spec)
if( !is.null(opt$help)||is.null(opt$input)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}
#default
if (is.null(opt$test)) {opt$test=0}
#Disable scientific notation
options("scipen"=100, "digits"=10)              
args <- commandArgs(T)
output <- opt$output
##RPKM
Data <- read.table(file=opt$input,header=F)           
Data <- as.matrix(Data)

#Prepare for gene expression quantification
search_gene_Data <- Data                          
search_gene_Data[1,1] <- "ID"
colnames(search_gene_Data) <- search_gene_Data[1,]
search_gene_Data <- search_gene_Data[-1,]

colnames(Data) <- Data[1,]
gene_col <- colnames(Data)
Data <- Data[-1,]
gene_id <- Data[,1]
rownames(Data) <- gene_id
Data <- na.omit(Data)
Data <- Data[,-1]
Data <- t(Data)
NewData <- matrix(nrow=nrow(Data),ncol=ncol(Data))
Data <- apply(Data,2,as.numeric)
Data <- log2(Data+1)
for(i in 1:ncol(Data)){
  if(mean(Data[,i],na.rm=TRUE)>=1)
    NewData[,i] <- Data[,i]
}
NewData <- t(NewData)
rownames(NewData) <- gene_id
colnames(NewData) <- gene_col[-1]
NewData <- na.omit(NewData)
#Filter expression levels
write.table(NewData,file = paste0(output,"/guolv.txt"),row.names=T,col.names=T,quote=F)    
ge <- rownames(NewData)
XgData <- matrix(nrow=nrow(NewData),ncol=nrow(NewData))
NewData <- t(NewData)
XgData <- abs(cor(NewData))
rownames(XgData) <- ge
#Calculate correlation coefficient
write.table(XgData,file = paste0(output,"/XgData.txt"),row.names=T,col.names=F,quote=F)  
###！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！###
Data <- read.table(file = paste0(output,"/XgData.txt"),header=F)
Data <- as.data.frame(Data)
NewData <- Data[,2:ncol(Data)]
rownames(NewData) = Data[,1]
NewData = apply(NewData,2,as.numeric)
NewData <- 1.-NewData
pp_thre = opt$td
thre = 1.-as.numeric(pp_thre)
dis <- as.dist(NewData)
hc <- hclust(dis,method="average")
cc <- cutree(hc,h=thre)
cc <- as.matrix(cc)
rownames(cc) = Data[,1]
#cluster
write.table(cc,file = paste0(output,"/class.txt"),row.names=T,col.names=F,quote=F)     
###！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！###
class_total <- read.table(file = paste0(output,"/class.txt"),header=F)
#Given a gene name
search_name <- opt$gene                                                                                          
for(i in 1:nrow(class_total)){
  if(class_total[i,1] == search_name)
  search_class <- class_total[i,2]
}
search_gene <- class_total
for(i in 1:nrow(search_gene)){
  if(search_gene[i,2] != search_class)
  search_gene[i,1] <- NA
}
search_gene <- na.omit(search_gene)
colnames(search_gene) <- c("ID","class")
search_gene_id <- search_gene[,1]
#Extract homologous genes
write.table(search_gene_id,file = paste0(output,"/search_gene.txt"),row.names=F,col.names=F,quote=F)      
###！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！###
if(opt$expression){
  gene_expression <- merge(search_gene,search_gene_Data,by.x="ID",by.y="ID")
  gene_expression <- gene_expression[,-2]
  write.table(gene_expression,file = paste0(output,"/search_gene_expression.txt"),row.names=F,col.names=F,quote=F)
}