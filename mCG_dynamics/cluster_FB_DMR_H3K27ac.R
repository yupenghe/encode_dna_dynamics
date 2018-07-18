tissues = c("FB")
stages <- c("E10_5","E11_5","E12_5",
            "E13_5","E14_5",
            "E15_5","E16_5","P0",
            "AD_1wk","AD_2wk","AD_4wk","AD_6wk")

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


        ## H3K27ac
        data <- read.table(paste0("results/",t,"DMR/",t,"DMR_",tissue,".tsv.H3K27ac"),
                           header=T,stringsAsFactors=F,sep="\t",comment = "$");
        y = data.matrix(data[,grep(tissue,colnames(data))])
        y[y < 0] = 0
        y <- apply(y,1,
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
        y = t(y)
        
        ## Exclude NA
        z = cbind(x,y)
        ind_nona = ind_nona = (apply(z,1,function(x) {return(sum(is.na(x)))}) == 0)
        z = z[ind_nona,]
        dmr = data[ind_nona,1:4]

        ## Clustering
        k = kmeans(z,centers=8,iter=100)
        m = aggregate(z[,1],list(k$cluster),mean)
        s = sort(m[,ncol(m)],decreasing=F,index.return=T)
        ind = NULL

        ## Save DMRs in each cluster
        write.table(            
            cbind(dmr,k$cluster),
            paste0(t,"DMR.bed"),
            quote=F,row.names=F,col.names=F,sep="\t")
    }
}
