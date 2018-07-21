library(amap)
library(gplots)
library(RColorBrewer)
library(fastcluster)
#source("/gale/netapp/home/yupeng/ENCODE/integrative/figures/data/utilities.R")

hclust_method = "ward"
dist_method= "euclidean"
#col_fun = colorRampPalette(c("dark blue",'light blue','yellow','red'))

pdf("heatmap_ENCODE_mCH_all_query.pdf",onefile=T,width=12)

mCH <- read.table("mCHdomain_merged.bed.profile",
                  header=T,stringsAsFactors=F,sep="\t");
mCH = mCH[,-c(grep("AD_HT",colnames(mCH)),
    grep("AD_KD",colnames(mCH)))]

num_up_bin = 20
num_domain_bin = 10
num_down_bin = 20
num_bin = num_up_bin + num_domain_bin + num_down_bin
        
## Heatmap
y = NULL
#for(ind in 4:ncol(mCH)){
for(ind in unique(c(grep("FB",colnames(mCH)),
                    grep("MB",colnames(mCH)),
                    grep("HB",colnames(mCH)),
                    grep("HT",colnames(mCH)),
                    grep("IT",colnames(mCH)),
                    grep("ST",colnames(mCH)),
                    grep("KD",colnames(mCH))#,
                    #grep("CF",colnames(mCH)),
                    #grep("LM",colnames(mCH)),
                    #grep("LG",colnames(mCH))
                    )
                  )){
    y = cbind(y,
        t(sapply(strsplit(mCH[,ind],","),as.numeric)))
    ##data.matrix(unlist(strsplit(mCH[,ind],","))))
}

y[is.na(y)] = 0

print(dim(y))

col_tissues <- c("firebrick3","lightsalmon4","coral","darkgoldenrod1",
                 "dodgerblue4","deepskyblue","cadetblue3","cornflowerblue",#"darkslategray2",
                 "darkseagreen","darkolivegreen","forestgreen",
                 "darkkhaki"
                 );
col_tissues = sample(col_tissues)

## Normalized
z = y
x = NULL
for(ind in 1:(ncol(y) / num_bin)){
    for(jnd in 1:nrow(y)){
        
        z[jnd,1:num_bin + (ind-1) * num_bin] = (y[jnd,1:num_bin + (ind-1) * num_bin] + 0.001) / (mean(y[jnd,c(1:num_up_bin,(num_bin - num_down_bin + 1):num_bin) + (ind-1) * num_bin],na.rm=T) + 0.001)
    }
    #x = cbind(x, z[,21:30 + (ind-1) * 30])
}

z = log2(z)
#x = log2(x)

rownames(y) = paste(mCH[,1],mCH[,2],mCH[,3],sep="-")
rownames(z) = paste(mCH[,1],mCH[,2],mCH[,3],sep="-")

ind = NULL
ind_clst = NULL
for(c in 1:5){
    domain = read.table(paste0("mCH_domain_g",c,".bed"),header=F,stringsAsFactors=F)
    ind = c(ind,paste(domain[,1],domain[,2],domain[,3],sep="-"))
    ind_clst = c(ind_clst,rep(c,nrow(domain)))
}


par(oma=c(1,1,1,3))
col_fun <- colorRampPalette(
    c(
        rev(brewer.pal(8,"RdYlBu"))[1:4],
        "white",
        rev(brewer.pal(8,"RdYlBu"))[5:8]
        )
    )

heatmap.2(z[ind,],trace='n',Colv=F,Rowv=F,
          key.xlab="normalized mCH level",
          density.info = 'none',key.title=NA,
          breaks=seq(-0.5,0.5,length.out=51),
          labCol = NA,
          labRow = NA,
          col=col_fun(50),
          RowSideColors = col_tissues[ind_clst],

          colsep = seq(0,ncol(y),by=num_bin),
          rowsep = c(0,as.numeric(cumsum(table(ind_clst)))),
          sepwidth=c(0.005,0.01),
          
          distfun=function(x) Dist(x,method=dist_method),
          hclustfun=function(x) hclust(x,method=hclust_method),
          sepcolor = "black"
          )

par(oma=c(1,1,1,3))

col_fun <- colorRampPalette(
    c(
        ##rev(brewer.pal(8,"RdYlBu"))[1:4],
        rev(brewer.pal(8,"RdYlBu"))[1:2],
                                        #"white",
        rev(brewer.pal(8,"RdYlBu"))[5:8]
        )
    )

h <- heatmap.2(y[ind,],trace='n',Colv=F,Rowv=F,
               key.xlab="mCH level",
               density.info = 'none',key.title=NA,
               breaks=seq(0,0.005,length.out=51),
               RowSideColors = col_tissues[ind_clst],
               labCol = NA,
               labRow = NA,
               col=col_fun(50),

               colsep = seq(0,ncol(y),by=num_bin),
               rowsep = c(0,as.numeric(cumsum(table(ind_clst)))),
               sepwidth=c(0.005,0.01),
               
               distfun=function(x) Dist(x,method=dist_method),
               hclustfun=function(x) hclust(x,method=hclust_method),
               sepcolor = "black"
               )



dev.off()
