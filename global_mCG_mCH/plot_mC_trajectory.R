#!/bin/env Rscript
library(ggplot2)
library(RColorBrewer)

col_fun = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))

df <- read.table("dat_sample_table.tsv",
                 header=T,
                 stringsAsFactors=F)

## Replace AD with p8w if possible
stages <- c("E10_5","E11_5",
            "E12_5","E13_5",
            "E14_5","E15_5",
            "E16_5","P0","AD");


## Reorder stages
stage_index = rep(NA,nrow(df))
for(i in 1:length(stages)){
    stage_index[df$stage == stages[i]] = i
}
df$stage = factor(df$stage,levels = df$stage[order(stage_index)])


tissues <- c("FB","MB","HB","NT",
             "HT","CF","LM","KD",
             "LG","ST","IT","LV");

## Reorder by tissues
tissue_index = rep(NA,nrow(df))
for(i in 1:length(tissues)){
    tissue_index[df$tissue == tissues[i]] = i
}
df$tissue = factor(df$tissue,levels = df$tissue[order(tissue_index)])

col_tissues <- c("firebrick3","lightsalmon4","coral","darkgoldenrod1",
                 "dodgerblue4","deepskyblue","cadetblue3","cornflowerblue",#"darkslategray2",
                 "darkseagreen","darkolivegreen","forestgreen",
                 "darkkhaki"
                 );

## mCG level
g <- ggplot(df) + aes(stage,mCG,group=tissue,color=tissue) +
    geom_line(size=1.1) +
        coord_cartesian(ylim = c(0.55,0.85)) +
            scale_color_manual(values=col_tissues)
                                

ggsave("trajectory_mCG.pdf",width=4, height=7)

## mCH level
g <- ggplot(df) + aes(stage,mCH,group=tissue,color=tissue) +
    geom_line(size=1.1) +
        coord_cartesian(ylim = c(0,0.0025)) +
        #coord_cartesian(ylim = c(0,0.015)) +
        #ylim(0,0.015) + 
            scale_color_manual(values=col_tissues)
                                

ggsave("trajectory_mCH.pdf",width=4, height=7)

## mCH level with large range
g <- ggplot(df) + aes(stage,mCH,group=tissue,color=tissue) +
    geom_line(size=1.1) +
        coord_cartesian(ylim = c(0,0.015)) +
            scale_color_manual(values=col_tissues)
                                

ggsave("trajectory_mCH_AD.pdf",width=4, height=7)

