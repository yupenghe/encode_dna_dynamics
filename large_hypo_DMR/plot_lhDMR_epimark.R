tissues <- c("FB","MB","HB","NT",
             "HT","LM","CF",
             "KD","ST","IT","LG");

pdf("lhDMR_epimark.pdf",onefile=T)
layout(t(1:2))
for(mark in c("H3K27ac","H3K4me1")){
    y_all = list(c(),c())
    for(tissue in tissues){
        input.file = paste0("results/lhDMR/lhDMR_",tissue,".region_with_epimark.tsv")
        lhdmr = read.table(input.file,header=T,stringsAsFactors=F)
        lhdmr = lhdmr[,grep(tissue,colnames(lhdmr))]
        lhdmr = lhdmr[,grep(mark,colnames(lhdmr))]

        input.file = paste0("results/lhDMR/non_lhDMR_",tissue,".region_with_epimark.tsv")
        non.lhdmr = read.table(input.file,header=T,stringsAsFactors=F)
        non.lhdmr = non.lhdmr[,grep(tissue,colnames(non.lhdmr))]
        non.lhdmr = non.lhdmr[,grep(mark,colnames(non.lhdmr))]

        y = list()
        for(s in colnames(lhdmr)){
            y = append(y,list(lhdmr[,s]))
            y = append(y,list(non.lhdmr[,s]))
            y_all[[1]] = c(y_all[[1]],lhdmr[,s])
            y_all[[2]] = c(y_all[[2]],non.lhdmr[,s])
        }
    }

    boxplot(y_all,
            names = c("lhDMR","non lhDMR"),
            notch=T,
            outline=T,
            main=mark,pch=16,cex=0.5,
            col=c("#CEA1A1","#809E84")
            )

}
dev.off()
