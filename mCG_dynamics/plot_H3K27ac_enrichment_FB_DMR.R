library(RColorBrewer)

stages <- c("E10_5","E11_5","E12_5",
            "E13_5","E14_5",
            "E15_5","E16_5","P0",
            "AD_1wk","AD_2wk","AD_4wk","AD_6wk")

layout(1:2)
tissue = "FB"

data <- read.table(paste0("results/devDMR/devDMR_",tissue,".tsv.Meth"),
                   header=T,stringsAsFactors=F,sep="\t",comment = "$");
colnames(data) = gsub("methylation_level_","",colnames(data))
x = data.matrix(data[,grep(tissue,colnames(data))])
y = NULL
for(stage in stages){
    ind = grep(paste0(stage,"_",tissue), colnames(x))
    y = cbind(y,x[,ind])
}
x = y
colnames(x) = paste0(stages,"_",tissue)
rownames(x) = data[,4]

## H3K27ac
data <- read.table(paste0("results/devDMR/devDMR_",tissue,".tsv.H3K27ac"),
                   header=T,stringsAsFactors=F,sep="\t",comment = "$");
y = data.matrix(data[,grep(tissue,colnames(data))])
y[y < 0] = 0
rownames(y) = data[,4]


## Exclude NA
z = cbind(x,y)
ind_nona = ind_nona = (apply(z,1,function(x) {return(sum(is.na(x)))}) == 0)
z = z[ind_nona,]
dmr = data[ind_nona,1:4]

stages <- c("E11_5","E12_5",
            "E13_5","E14_5",
            "E15_5","E16_5","P0")
y = NULL
for(stage in stages){
    mcg = z[,paste0(stage,"_FB")]
    k27ac = z[,paste0("H3K27ac_",stage,"_FB")]

    ind = mcg <= 0.5

    y <- cbind(y, c(sum(ind & k27ac >=0 & k27ac <=2 ),
                    sum(ind & k27ac >2 & k27ac <=4 ),
                    sum(ind & k27ac >4 & k27ac <=6 ),
                    sum(ind & k27ac >6 )) / sum(ind) )

    y <- cbind(y, rep(NA,4))
    
    y <- cbind(y, c(sum(!ind & k27ac >=0 & k27ac <=2 ),
                    sum(!ind & k27ac >2 & k27ac <=4 ),
                    sum(!ind & k27ac >4 & k27ac <=6 ),
                    sum(!ind & k27ac >6 )) / sum(!ind) )
    y <- cbind(y, rep(NA,4))
    y <- cbind(y, rep(NA,4))

}



pdf("Figure_2f.pdf",onefile=T,width=7)
col_fun <- colorRampPalette(rev(brewer.pal(8,"RdYlBu")))
barplot(y,col = col_fun(4))

legend(1,1,c("> 6",
             "4 ~ 6",
             "2 ~ 4",
             "0 ~ 2"),
       fill = col_fun(4))

dev.off()
