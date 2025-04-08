# Comprehensive Multi-Omic Insights into High-Grade Glioma in Children, Adolescents and Young Adults 
Nicole L. Tignor*, Mateusz Koptyra*2, Shrabanti Chowdhury*1, Francesca Petralia*1, Marina A. Gritsenko*3, Weiping Ma*1, Yuankun Zhu**2, Giacomo B. Marino**1, Xiaoyu Song**4, Simona Migliozzi**5, Jeffrey R. Whiteaker**6, Dmitry Rykunov**1, Yingwei Hu**7, Noshad Hosseini8, Jo Lynne Rokita9,10, Felipe da Veiga Leprevost11, Komal S. Rathi2, Lijun Chien7, Yi-Ting Wang12, Karl K. Weitz12, Rosalie K. Chu3, Ronald J. Moore12, Azra Krek1, Xuran Wang1, Eden Z. Deng1, Chia-Feng Tsai12, Tyler Sagendorf12, Vladislav A. Petyuk12, Tujin Shi12, Thomas L. Fillmore3, Rui Zhao12, Matthew E. Monroe12, Marius V. Dannappel13, Paul Daniel13, Lei Zhao6, Richard G. Ivey6, Uliana J. Voytovich6, Tomer M. Yaron-Barir14,15, Emily M. Huntsman14, Jared L. Johnson16, Nakib Abedin1, Yan-chak Li1, Mariarita Santi17, Demetri Dupal2, Jena Lilly2, Adam Kraya2, Joseph M. Dybas2, Bo Zhang2, Chuwei Zhong17, Miguel A. Brown2, Saksham Phul17, Eric Wafula18, Alvin Farrel18, Zhuangzhuang Geng2, Ryan J. Corbett2, Ammar S. Naqvi2, Daniel P. Miller2, Jennifer Mason2, Tatiana S. Patton2, Stephanie McGrory2, Shannon Robins2, Allison Heath2, Catherine Sullivan2, Noel Coleman2, Allison Morgan2, Luciano Garofano5, Boris Reva1, Eric E. Schadt1, Richard D. Smith12, Mehdi Mesri19, Ana I. Robles19, Lewis C. Cantley16, Li Ding20, Karin D. Rodland21, Bing Zhang22, Alexey Nesvizhskii11,8, Antonio Iavarone5, Marcin Cieslik11,8, Joseph Ippolito23, Phillip B. Storm‡2, Josh Rubin‡24, Ron Firestein‡13, Avi Ma'ayan‡1, Hui Zhang‡7,25, Amanda Paulovich‡6, Tao Liu‡12, Adam Resnick††2, Brian Rood††9, Pei Wang††1,26, Philadelphia Coalition for a Cure, Children's Brain Tumor Network, Clinical Proteomic Tumor Analysis Consortium
* Co-first authors (equal contributors).  
** Co-second authors (equal contributors).
‡ Co-senior authors
†† Co-corresponding authors 
26lead contact: pei.wang@mssm.edu


Module authors: [Komal S. Rathi](https://github.com/komalsrathi/),
                [Zhuangzhuang Geng](https://github.com/zzgeng),
                [Jo Lynne Rokita](https://github.com/jharenza),
		[Joseph M. Dybas](https://github.com/JosephDybas)

## Clone the repository

To clone the repository, run:
```
git clone git@github.com:d3b-center/hope-cohort-analysis.git
```

## Download the data

To download the current data release:
```
bash download-data.sh
```

## Run the script on docker

To pull the docker image, run the command line:
```
docker pull pgc-images.sbgenomics.com/d3b-bixu/d3b_hope_analysis:latest
```

To start the docker container, run the command line below. For mac M1 user, add `--platform=linux/arm64`.
```
docker run -d -e PASSWORD=pass -p 8787:8787 --name <CONTAINER_NAME> -v $PWD:/home/rstudio/hope-cohort-analysis pgc-images.sbgenomics.com/d3b-bixu/d3b_hope_analysis:latest
```

To use docker in command line:
```
docker exec -ti <CONTAINER_NAME> bash
```

## Navigate to the repository root
```
cd /home/rstudio/hope-cohort-analysis
```

## Modules

```
analyses
├── alt-analysis
├── data-availability
├── master-annotation
├── merge-files
├── tp53_nf1_score
├── molecular-subtyping-HGG
├── msi-sensor-analysis
├── oncoplots
└── survival-analysis 
```

1) `data-availability`: This module has scripts to create data availability plots.
2) `merge-files`: This module has scripts to merge files obtained from cavatica i.e. RSEM gene expression, Consensus MAF, ControlFREEC, Fusions which are then filtered and annotated. 
3) `master-annotation`: This module combines various sources of information from the HOPE group into one single tsv file for downstream analyses.
3) `msi-sensor-analysis`: Downstream analyses with MSISensor pro outputs.
4) `oncoplots`: This module has scripts to create oncoplots and cascade plots. Reference files and genelists were obtained from [PNOC003](https://github.com/d3b-center/d3b-pnoc003-HGG-DMG-omics/tree/master/analyses/Oncoplot)
5) `survival-analysis`: This module has scripts to do survival analysis with ALT status and molecular subtypes.
6) `alt-analysis`: Downstream analyses with ALT status. 
7) `tmb-calculation`: Adapted from [OpenPedCan-anaysis](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/tmb-calculation)
