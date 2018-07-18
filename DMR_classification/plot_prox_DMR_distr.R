x = c(
    ##Proximal CGI-promoter DMR
    46692,
    ##Proximal nonCGI-promoter DMR
    90831,
    ##Proximal CGI DMR
    1710,
    ##Proximal CGI-shore DMR
    13786)
labels = c("CGi promoter",
    "non-CGi promoter",
    "CGi",
    "CGi shore")
names(x) = paste0(labels,"\n(",
         round(x/sum(x)*100,digits=1),
         "%, n = ",x,")")

pdf("Extended_Data_Fig_2b.pdf")
pie(x)
