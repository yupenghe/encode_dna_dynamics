# Spatiotemporal DNA Methylome Dynamics of the Developing Mammalian Fetus

* Code collection and depositing on this site is in progress.
* All customs script will be deposited within 1 week.

### Data collection
Scripts to download, from ENCODE DCC, the ChIP-seq and ATAC-seq processed bam files of E10.5-P0 samples.
```bahs
cd data_collection/
python get_files_from_DCC_ATAC.py
python get_files_from_DCC_ChIP.py
```
Running these two scripts will take ~1T space and several hours. 

### Global mCG and mCH levels
Plot the global mCG and mCH level of samples
```bash
Rscript plot_mC_trajectory.R
```
