library(amap)
library(gplots)
library(RColorBrewer)
library(fastcluster)
#source("/gale/netapp/home/yupeng/ENCODE/integrative/figures/data/utilities.R")

hclust_method = "ward"
dist_method= "euclidean"
#col_fun = colorRampPalette(c("dark blue",'light blue','yellow','red'))
col_fun <- colorRampPalette(
    c(
        rev(brewer.pal(8,"RdYlBu"))[1:4],
        "white",
        rev(brewer.pal(8,"RdYlBu"))[5:8]
        )
    )

pdf("heatmap_ENCODE_mCH_all_tp.pdf",onefile=T)

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
for(ind in unique(c(#grep("FB",colnames(mCH)),
                    grep("MB",colnames(mCH)),
                    grep("HB",colnames(mCH)),
                    grep("HT",colnames(mCH)),
    grep("IT",colnames(mCH)),
                    grep("ST",colnames(mCH)),
    grep("KD",colnames(mCH))
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

k = kmeans(z,centers=5,iter.max = 100)
s = sort(k$cluster,index.return=T)

for(c in unique(k$cluster)){
    write.table(mCH[k$cluster == c,1:3],
                paste0("mCH_domain_g",c,".bed"),
                quote=F,row.names=F,col.names=F,
                sep="\t")
}


par(oma=c(1,1,1,3))
heatmap.2(z[s$ix,],trace='n',Colv=F,Rowv=F,
          key.xlab="mCH/CH",
          density.info = 'none',key.title=NA,
          breaks=seq(-0.5,0.5,length.out=51),
          labCol = NA,
          labRow = NA,
          col=col_fun(50),
          RowSideColors = col_tissues[s$x],

          colsep = seq(0,ncol(y),by=num_bin),
          rowsep = c(1,as.numeric(cumsum(table(s$x)))),
          sepwidth=c(0.01,0.05),
          
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

h <- heatmap.2(y[s$ix,],trace='n',Colv=F,Rowv=F,
               key.xlab="mCH/CH",
               density.info = 'none',key.title=NA,
               breaks=seq(0,0.005,length.out=51),
               RowSideColors = col_tissues[s$x],
               labCol = NA,
               labRow = NA,
               col=col_fun(50),

               colsep = seq(0,ncol(y),by=num_bin),
               rowsep = c(1,as.numeric(cumsum(table(s$x)))),
               sepwidth=c(0.01,0.05),
               
               distfun=function(x) Dist(x,method=dist_method),
               hclustfun=function(x) hclust(x,method=hclust_method),
               sepcolor = "black"
               )



dev.off()
