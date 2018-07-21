#library("changepoint")

NCR <- read.table("../data/stats/all_NCR.tsv",
                  header=T,stringsAsFactors=F)
ncr = NCR[,2]
names(ncr) = NCR[,1]

nsite <- read.table("data/mm10_bin5kb_with_CH_nsite.tsv",header=T,
                    stringsAsFactors=F)
colnames(nsite) = gsub("nsite_","",colnames(nsite))


mCH <- read.table("data/mm10_bin5kb_with_CH_level.tsv",header=T,
                  stringsAsFactors=F);

mCH = mCH[,grep("AD",colnames(mCH),invert=T)]
colnames(mCH) = gsub("methylation_level_","",colnames(mCH))
dmr = mCH[,1:3]


for(s in colnames(mCH)[4:ncol(mCH)]){

    mCH[,s] = mCH[,s] - ncr[s]
    mCH[mCH[,s] < 0 & !is.na(mCH[,s]) ,s] = 0
}

## Filter out regions that are susceptible to mapping error

if( !(is.na(file.info("tmp/excl_mCHdomain.bed")$size) |
          file.info("tmp/excl_mCHdomain.bed")$size == 0)
   ){
    
    excl_domain <- read.table("tmp/excl_mCHdomain.bed",header=F,
                              stringsAsFactors=F)
    for(i in 1:nrow(excl_domain)){
        ind = mCH$chr == excl_domain[i,1] &
            mCH$start >= excl_domain[i,2] &
                mCH$end <= excl_domain[i,3]
        mCH[ind,4:ncol(mCH)] = NA
    }
}



## Begin
chromosomes = 1:19

suppressPackageStartupMessages(library("doParallel",,quietly=T,verbose=F))   
suppressPackageStartupMessages(library("foreach",quietly=T,verbose=F))

cl <- makeCluster(16)
registerDoParallel(cl)


args <- commandArgs(TRUE)

if(length(args) != 0 & args[1] == "Control"){
    sample_list <- c("E10_5_FB_1","E10_5_FB_2",
                     "E10_5_MB_1","E10_5_MB_2",
                     "E10_5_HB_1","E10_5_HB_2",
                     "E10_5_CF_1","E10_5_CF_2",
                     "E10_5_LM_1","E10_5_LM_2",
                     ##"E10_5_HT_1","E10_5_HT_2",
                     ##E11.5
                     "E11_5_CF_1","E11_5_CF_2",
                     "E11_5_FB_1","E11_5_FB_2",
                     "E11_5_HB_1","E11_5_HB_2",
                     ##"E11_5_HT_1","E11_5_HT_2",
                     "E11_5_LM_1","E11_5_LM_2",
                     "E11_5_MB_1","E11_5_MB_2",
                     "E11_5_NT_1","E11_5_NT_2",
                     "E11_5_LV_1","E11_5_LV_2");
}else{
    sample_list = colnames(mCH)[4:ncol(mCH)]
}

foreach(s = sample_list) %dopar% {
#for(s in colnames(mCH)[4:ncol(mCH)]) {
    print(s)
    library("changepoint")
    out = NULL
    for(chrom in chromosomes){
        print(chrom)
        ind = mCH[,1] == chrom & (!is.na(mCH[,s])) &
            nsite[,"mm10"] >= 500 &
                (nsite[,s]/(nsite[,"mm10"]) >= 0.5)
    
        domain = cpt.mean(mCH[ind,s]*1000,method="PELT",pen.value=0.01,penalty="Asymptotic",minseglen=2)
        ##domain = cpt.mean(mCH[ind,s],method="PELT",pen.value=0.05,penalty="Asymptotic")

        chgpt = which(ind)[cpts(domain)]

        if(length(chgpt) > 0) {
            coor = NULL
            M = NULL
            chgpt = c(1,chgpt)
            if(max(chgpt) < length(ind)) { chgpt = c(chgpt,length(ind)) }
            
            for(i in 2:length(chgpt)){
                M = c(M,mean(mCH[(chgpt[i-1]+1):chgpt[i],s],na.rm=T))
                coor = rbind(coor,t(c(chrom,mCH[chgpt[i-1],3],mCH[chgpt[i],3])))
            }
            out = rbind(out,cbind(coor,M))
        }
    }
    write.table(out,
                paste0("tmp/mCH_segment_",s,".bed"),
                quote=F,row.names=F,col.names=F,
                sep="\t")
}
