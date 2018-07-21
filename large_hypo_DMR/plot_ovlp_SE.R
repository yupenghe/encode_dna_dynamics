
se = read.table("superenhancer/ovlp_SE.tsv",stringsAsFactors=F)

pdf("bar_ovlp_gene.pdf",height=5)
y = data.matrix(se[,2:3])
y[,1] = y[,1] - y[,2]
y = y[,c(2,1)]
y = t(y)
b = barplot(y,names=se[,1],
    ylab = "Number of large hypo CG-DMRs",
    col = c("dark grey","light grey")
            )

print(y)
perc = paste0(round(y[1,]/(colSums(y))*100,digits=2),"%")
print(perc)
text(b,y[1,]+10,perc,
     xpd=T)

dev.off()
