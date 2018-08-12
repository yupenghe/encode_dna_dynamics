# Spatiotemporal DNA Methylome Dynamics of the Developing Mammalian Fetus

## Dependencies
* python2/3
* R
* (R package) riverplot
* (R package) ggplot2
* (R package) RColorBrewer
* (R package) amap
* (R package) gplots
* (R package) fastcluster
* (R package) doParallel
* (R package) foreach
* (R package) changepoint

## Scripts
### Data collection
Scripts to download, from ENCODE DCC, the ChIP-seq and ATAC-seq processed bam files of E10.5-P0 samples.
```bahs
cd data_collection/
python get_files_from_DCC_ATAC.py
python get_files_from_DCC_ChIP.py
```
Running these two scripts will take ~1T space and several hours. 

### Global mCG and mCH levels
Plot the global mCG and mCH level of samples (Figure 1b and Figure 4a).
```bash
cd global_mCG_mCH/
Rscript plot_mC_trajectory.R
```
### CG-DMR classification
Riverplot showing CG-DMR classification (Figure 1e).
```bash
Rscript DMR_classification/plot_DMR_classification.R
```

Pie chart of proximal CG-DMRs (Extended Data Figure 2b)
```bash
Rscript DMR_classification/plot_prox_DMR_distr.R
```

### mCG dynamics
mCG dynamics in tissue-specific CG-DMRs (Figure 2a-d).
```bash
Rscript mCG_dynamics/plot_barchart_merge.R
Rscript mCG_dynamics/plot_Heatmap_merge.R
```

`mCG_dynamics/cluster_FB_DMR_H3K27ac.R` is used to cluster forebrain-specific CG-DMRs
with distinct mCG and H3K27ac dynamics during development into 8 groups. Two scripts
below plot the mCG and H3K27ac profile of different classes of forebrain-specific CG-DMRs, and the correlation between mCG and H3K27ac in these CG-DMRs (Figure 2e-f).
```bash
Rscript mCG_dynamics/plot_Heatmap_H3K27ac_FB.R
Rscript mCG_dynamics/plot_H3K27ac_enrichment_FB_DMR.R
```

### large hypo- CG-DMR
`large_hypo_DMR/plot_lhDMR_epimark.R` and `large_hypo_DMR/plot_ovlp_SE.R` are used
to plot the intensity of epigenetic modification in large hypo- CG-DMRs, and the overlap
between large hypo- CG-DMRs and super-enhancers (Figure 3b, c). 

### mCH domain calling
`mCH_domain_calling/get_mCHdomain.pl` and `mCH_domain_calling/call_changepoint.R` are scripts used to call mCH domains.

### Clustering mCH domains
Scripts in `mCH_domain_clustering/` are used to cluster mCH domains based on their mCH dynamics across tissues as well as visualization of the mCH dynamics of clustered mCH domain (Figure 4c and Extended Data Figure 7d).

