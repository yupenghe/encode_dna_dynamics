# Spatiotemporal DNA Methylome Dynamics of the Developing Mammalian Fetus

* Code collection and depositing on this site is in progress.
* All customs script will be deposited within 1 week.

## Dependencies
* python2/3
* R
* (R package) riverplot
* (R package) ggplot2
* (R package) RColorBrewer
* (R package) amap
* (R package) gplots
* (R package) fastcluster

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
