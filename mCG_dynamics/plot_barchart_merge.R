stages <- c("E10_5","E11_5","E12_5",
            "E13_5","E14_5","E15_5",
            "E16_5","P0","AD")

tissues <- c("FB","MB","HB","NT",
             "HT","LM","CF",
             "KD","ST","IT","LG");

pdf("Figure_2b.pdf",onefile=T,width=3)
z_all = NULL
layout(matrix(1:22,byrow=F,nrow=11))
for(tissue in tissues){
    for(t in c("dev")){
        data <- read.table(paste0("results/",t,"DMR/",t,"DMR_",tissue,".tsv.Meth"),
                           header=T,stringsAsFactors=F,sep="\t",comment = "$");
        data[,1] = paste0("chr",data[,1])
        colnames(data) = gsub("methylation_level_","",colnames(data))

        x = data.matrix(data[,grep(tissue,colnames(data))])

        y = NULL
        ind = NULL
        for(i in 2:length(stages)){
            stage = stages[i]
            if(tissue == "FB" & stage == "AD"){
                stage = "AD_6wk"
            }
            if(tissue == "HT" | tissue == "KD" | tissue == "LG" ){
                if(stage == "AD"){
                    stage = "p8w"
                }
            }

            if(length(grep(paste0(stage,"_",tissue),colnames(x)))==0 |
               length(grep(paste0(stages[i-1],"_",tissue),colnames(x)))==0){
                y = cbind(y,rep(NA,nrow(x)))
                ind = c(ind,F)
            }else{
                meth_diff = x[,grep(paste0(stage,"_",tissue),colnames(x))] -
                    x[,grep(paste0(stages[i-1],"_",tissue),colnames(x))]
                y = cbind(y,meth_diff)
                ind = c(ind,T)
                ## Write to files
                write.table(            
                    na.omit(data[meth_diff >= 0.1 ,1:4]),
                    paste0("results/stagewise/devDMR_hyper_",tissue,".",
                           stages[i-1],".",stage,
                           ".bed"),
                    quote=F,row.names=F,col.names=F,sep="\t")
                write.table(            
                    na.omit(data[meth_diff <= -0.1 ,1:4]),
                    paste0("results/stagewise/devDMR_hypo_",tissue,".",
                           stages[i-1],".",stage,
                           ".bed"),
                    quote=F,row.names=F,col.names=F,sep="\t")
                
            }
        }
        
        par(mar=c(1,3,1,3))

        z <- rbind(apply(y,2,function(x) {return( sum(x <= -0.1,na.rm=T))}),
                   apply(y,2,function(x) {return( sum(x >= 0.1,na.rm=T))}))
        z[,!ind] = NA

        if(t == "dev"){
            cat(tissue,nrow(data),"loss",z[1,],"\n",sep="\t")
            cat(tissue,nrow(data),"gain",z[2,],"\n",sep="\t")
        }

        
        z_all = rbind(z_all,unname(z)/nrow(x))
        b <- barplot(unname(z),beside=T,col=c("lightblue","dark red"),border=NA,
                     ylim = c(floor(min(z,na.rm=T)/10000)*10000,ceiling(max(z,na.rm=T)/10000)*10000)
                     )
        
        z <- rbind(apply(y,2,function(x) {return( sum(x <= -0.2,na.rm=T))}),
                   apply(y,2,function(x) {return( sum(x >= 0.2,na.rm=T))}))
    }
}
dev.off()

pdf("Figure_2c_d.pdf")
layout(matrix(1:4,byrow=T,nrow=2))
par(mar=c(1,3,1,3))
for(t in 1:1){
    z_hypo = z_all[seq(1+(t-1)*2,nrow(z_all),by=2),]
    z_hyper = z_all[seq(2+(t-1)*2,nrow(z_all),by=2),]
    
    plot(c(1,ncol(z_all)),c(0,1),
         ylim=c(0,1),
         type='n')
    for(i in 1:length(tissues)){
        lines(1:ncol(z_all),z_hypo[i,],col="grey")
    }
    lines(1:ncol(z_all),colMeans(z_hypo,na.rm=T),col="blue")
    
    plot(c(1,ncol(z_all)),c(0,1),
         ylim=c(0,1),
         type='n')
    for(i in 1:length(tissues)){
        lines(1:ncol(z_all),z_hyper[i,],col="grey")
    }
    lines(1:ncol(z_all),colMeans(z_hyper,na.rm=T),col="red")
}
dev.off()
