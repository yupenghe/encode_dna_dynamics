library(riverplot)

pdf("Figure_1e.pdf")
nodes = data.frame(
    c("Total","Distal","Proximal",
      "fe","feflk","Remain",
      "cls",'uncls',"te","un"),
    c(1,2,2,
      3,3,3,
      4,4,5,5),
    rep(NA,10)
    )
names(nodes) = c("ID","x","labels")

edges = data.frame(
    c("Total","Total",
      rep("Distal",3),
      rep("Remain",2),
      rep("uncls",2)),
    c("Distal","Proximal",
      "fe","feflk","Remain",
      "cls",'uncls',"te","un"),
    c(1655791,153019,
      397320,212620,1045851,
      159347,886504,449623,436881)
    )
names(edges) = c("N1","N2","Value")
plot(makeRiver(nodes,edges),
     plot_area = 0.8, yscale=0.06)
dev.off()
