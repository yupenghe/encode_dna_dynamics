normalize <- function(y){
    z <- apply(y,1,
               function(x) {
                   x_std = sd(x)
                   if(x_std == 0){
                       return(rep(0,length(x)))
                   }else{
                       #return((x - mean(x))/x_std)
                       return( x / max(x))
                   }
               }
               )
    return(t(z))
}


library(amap)
library(gplots)
library(RColorBrewer)
library(fastcluster)
source("utilities.R")

hclust_method = "ward"
dist_method= "euclidean"
col_fun = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))

tissues = c("FB")
stages <- c("E10_5","E11_5","E12_5",
            "E13_5","E14_5",
            "E15_5","E16_5","P0",
            "AD_1wk","AD_2wk","AD_4wk","AD_6wk")

pdf("Figure_2e.pdf",onefile=T,width=7)
#layout(matrix(1:4,byrow=T,nrow=2))
for(tissue in tissues){
    for(t in c("dev")){
        data <- read.table(paste0("results/",t,"DMR/",t,"DMR_",tissue,".tsv.Meth"),
                           header=T,stringsAsFactors=F,sep="\t",comment = "$");
        colnames(data) = gsub("methylation_level_","",colnames(data))
        x = data.matrix(data[,grep(tissue,colnames(data))])
        y = NULL
        for(stage in stages){
            ind = grep(paste0(stage,"_",tissue), colnames(x))
            y = cbind(y,x[,ind])
        }
        x = y
        rownames(x) = data[,4]

        ## H3K27ac
        data <- read.table(paste0("results/",t,"DMR/",t,"DMR_",tissue,".tsv.H3K27ac"),
                           header=T,stringsAsFactors=F,sep="\t",comment = "$");
        y = data.matrix(data[,grep(tissue,colnames(data))])
        y[y < 0] = 0
        #y = normalize(y)
        rownames(y) = data[,4]

        ## RNA
        data <- read.table(paste0("results/",t,"DMR/",t,"DMR_",tissue,".tsv.rna"),
                           header=T,stringsAsFactors=F,sep="\t",comment = "$");
        tpm = data.matrix(data[, grep(tissue,colnames(data))])
        tpm = log10(tpm+1)
        #tpm = normalize(tpm)
        rownames(tpm) = data[,4]
        
        par(mar=c(2,0.5,2,0.5))

        ## Exclude NA
        z = cbind(x,y)
        ind_nona = (apply(z,1,function(x) {return(sum(is.na(x)))}) == 0)
        z = z[ind_nona,]
        dmr = data[ind_nona,1:4]
        tpm = tpm[ind_nona,]
        tpm_dist = abs(data$dist[ind_nona])
        ## Clustering

        ## read the DMRs in different clusters
        clst_dmr = read.table(paste0("DMR_FB/",t,"DMR.bed"),header=F,stringsAsFactors=F)
        rownames(clst_dmr) = clst_dmr[,4]
        clst_dmr = clst_dmr[rownames(z),5]
        m = aggregate(z[,1],list(clst_dmr),mean)
        s = sort(m[,ncol(m)],decreasing=F,index.return=T)
        print(s$ix)
        ind = NULL
        num_clst = max(clst_dmr)        
        ## Plot profiles
        layout(1:num_clst)
        for(c in 1:num_clst){
            ind = c(ind,which(clst_dmr == c))
        }        
        for(c in 1:num_clst){
            plot(c(1,14),c(0,1),type='n',
                 main=paste0("cluster ",c," - ",sum(clst_dmr == c)),
                 xaxt='n',xaxs='i'
                 )
            axis(1,1:14,NA)
            text(1:14,-0.1,
                 c("E10_5","E11_5","E12_5",
                   "E13_5","E14_5",
                   "E15_5","E16_5","P0",
                   "P1w","P2w","P3w","P4w","P6w","P8w"),
                 xpd=T,srt=270)
            
            lines(c(1:8,c(1,2,4,5)+8),colMeans(z[clst_dmr == c,1:ncol(x)]),col="blue")
            lines(c(2:8,c(1,3,6)+8),colMeans(z[clst_dmr == c,1:ncol(y)+ncol(x)]),col="red")
        }
        z = z[ind,]

        ## Plot Heatmaps
        layout(t(1:3))
    
        ## Methylation
        image(t(z[,1:ncol(x)]),col=col_fun(50),breaks=seq(0,1,length.out=51),
              axes=F)

        ## H3K27ac
        image(t(z[,ncol(x)+1:ncol(y)]),col=col_fun(50),
              breaks=c(seq(0,6,length.out=50),1e6),
              axes=F)
        
    }
}
dev.off()
