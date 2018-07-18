library(amap)
library(gplots)
library(RColorBrewer)
library(fastcluster)
source("utilities.R",local=TRUE)

hclust_method = "ward"
dist_method= "euclidean"
col_fun = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))

stages <- c("E10_5","E11_5","E12_5",
            "E13_5","E14_5",
            "E15_5","E16_5","P0","AD")

tissues <- c("FB","MB","HB","NT",
             "HT","LM","CF",
             "KD","ST","IT","LG");

pdf("Figure_2a",onefile=T,width=1.5)
layout(matrix(1:11,byrow=T,nrow=11))
for(tissue in tissues){    
    data <- read.table(paste0("results/devDMR/devDMR_",tissue,".tsv.Meth"),
                       header=T,stringsAsFactors=F,sep="\t",comment = "$");
    colnames(data) = gsub("methylation_level_","",colnames(data))
    x = data.matrix(data[,grep(tissue,colnames(data))])
    x = smooth_matrix(x,50,50)
    
    ##Sort
    k = kmeans(x,centers=10,iter=100)
    m = aggregate(x[,1],list(k$cluster),mean)
    s = sort(m[,ncol(m)],decreasing=F,index.return=T)
    ind = NULL
    for(c in s$ix){
        ind = c(ind,which(k$cluster == c))
    }
    x = x[ind,]
    y = NULL
    for(stage in stages){
        ind = grep(paste0(stage,"_",tissue), colnames(x))
        if(tissue == "FB" & stage == "AD"){
            ind = grep("AD_6wk_FB", colnames(x))
        }
        if(tissue == "HT" | tissue == "KD" | tissue == "LG" ){
            if(stage == "AD"){
                stage = "p8w"
            }
        }
            
        if(length(ind) == 0){
            y = cbind(y,rep(NA,nrow(x)))
        }else if (length(ind) == 1){
            y = cbind(y,x[,ind])
        }else{
            y = cbind(y,rowMeans(x[,ind]))
        }
        
        next;
        
        if(!(paste0(stage,"_",tissue) %in% colnames(x))){
            y = cbind(y,rep(NA,nrow(x)))
        }else{
            y = cbind(y,x[,paste0(stage,"_",tissue)])
        }
    }
    colnames(y) = paste0(stages,"_",tissue)
    
    par(mar=c(0.5,0.5,0.5,0.5))
    
    image(t(y),col=col_fun(50),breaks=seq(0,1,length.out=51),
          axes=F)
}
dev.off()
