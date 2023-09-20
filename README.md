# HOPE cohort analysis

Module authors: [Komal S. Rathi](https://github.com/komalsrathi/),
                [Zhuangzhuang Geng](https://github.com/zzgeng),
                [Jo Lynne Rokita](https://github.com/jharenza)

## Clone the repository

To clone the repository, run:
```
git clone git@github.com:d3b-center/hope-cohort-analysis.git
```

## Download the data

To download the current data release:
```
bash download_data.sh
```

## Run the script on docker

To pull the docker image, run the command line:
```
docker pull pgc-images.sbgenomics.com/zhuangzhuanggeng/d3b_hope_analysis:latest
```

To run the docker, run the command line below. For mac M1 user, add `--platform=linux/arm64`.
```
docker run -d -e PASSWORD=pass -p 8787:8787 --name <CONTAINER_NAME> -v $PWD:/home/rstudio/HOPE pgc-images.sbgenomics.com/zhuangzhuanggeng/d3b_hope_analysis:latest
```

## Navigate to the repository root
```
cd /home/rstudio/d3b_hope_analysis
```


## Modules

```
analyses
├── alt-analysis
├── data-availability
├── master-annotation
├── merge-files
├── molecular-subtyping-HGG
├── msi-sensor-analysis
├── oncoplots
├── survival-analysis
└── tp53_nf1_score
```

1) `data-availability`: This module has scripts to create data availability plots.
2) `merge-files`: This module has scripts to merge files obtained from cavatica i.e. RSEM gene expression, Consensus MAF, ControlFREEC, Fusions which are then filtered and annotated. 
3) `master-annotation`: This module combines various sources of information from the HOPE group into one single tsv file for downstream analyses.
3) `msi-sensor-analysis`: Downstream analyses with MSISensor pro outputs.
4) `oncoplots`: This module has scripts to create oncoplots and cascade plots. Reference files and genelists were obtained from [PNOC003](https://github.com/d3b-center/d3b-pnoc003-HGG-DMG-omics/tree/master/analyses/Oncoplot)
5) `survival-analysis`: This module has scripts to do survival analysis with ALT status and molecular subtypes.
6) `alt-analysis`: Downstream analyses with ALT status. 
7) `tmb-calculation`: Adapted from [OpenPedCan-anaysis](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/tmb-calculation)
