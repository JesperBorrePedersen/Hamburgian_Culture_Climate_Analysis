## Research Compendium for 'Climate Niche Modeling Reveals the Fate of Pioneering Late Pleistocene Populations in Northern Europe'.

### Compendium DOI

[![DOI](https://zenodo.org/badge/286960061.svg)](https://zenodo.org/badge/latestdoi/286960061)

### Maintainer of the repository:

[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--3468--0986-green.svg)](https://orcid.org/0000-0002-3468-0986) Jesper Borre Pedersen (<jesper-borre.pedersen@ifu.uni-tuebingen.de>) 

### Description

This repository contains the data, code and supplementary information for our manuscript: _Climate niche modelling reveals the fate of pioneering Late Pleistocene populations in northern Europe._

### Content
1. [Authors](#Authors)
2. [Citation](#Citation)
3. [Acknowledgements](#Acknowledgements)
4. [Getting Started](#Getting-Started-with-the-Code)
5. [Software Requirements](#Software-Requirements)
6. [Folder Structure](#folder-structure)
7. [Model Evaluation and Animated Time-Series](#model-evaluation-and-animated-time-series)
8. [Landscape Metrics](#landscape-metrics)
8. [License](#License)

### Authors
Jesper Borre Pedersen, Jakob Johan Assmann, Signe Normand, Dirk Nikolaus Karger, Felix Riede.

### Citation

Pedersen, Jesper B., Assmann, Jakob J., Normand, Signe, Karger, Dirk N. & Felix Riede (2023): _Climate niche modelling reveals the fate of pioneering Late Pleistocene populations in northern Europe_.Current Anthropology 64(5): 599–608. https://doi.org/10.1086/726700

### Acknowledgements
(from the manuscript)

JFR’s and JJA contribution to this paper is part of CLIOARCH, an ERC Consolidator Grant project that has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (grant agreement No. 817564). DNK acknowledges funding from: The WSL internal grant exCHELSA and ClimEx, BioDivERSa (projects 20BD21_184131, 20BD21_193907), as well as the Swiss Data Science Projects: SPEEDMIND, COMECO, the WSL internal grant ClimEx. JBP would like to thank the University of Aarhus for granting a PhD fellowship.

### Getting Started with the Code
_Please note:_ This repository contains the code and the archaeological data required to replicate our analysis. However, the by-century temperature and preciptation rasters from the CHELSA Trace21k dataset are also needed. These data form the basis of the model fits and predictions. __The CHELSA Trace21k dataset is currenlty in review__ and not yet publicly available. We will update this section once the data has completed peer-review and has been released.  

The folder structure of this repository and the software requirements for the anlysis are outlined below. 

All code required for the analysis and reproduction of the figures can be found in the `scripts/analysis.R` script. The code is structured using RStudio's section delination syntax to aid navigation within the script. 

### Software Requirements
The analysis was developed and carried out using R version 3.6.0 and the following packages: sf 0.9-3, raster 3.1-5, tidyverse 1.3.0, cowplot 1.0.0, rasterVis 0.47, colorspace 1.4-1, dismo 1.3-3, magick 2.5.2 and landscapemetrics 1.5.1 .

### Folder Structure

```
├── data            Archaeological data       
├── figures         Figure and animation outputs
├── scripts         The analysis script(s)
└── tables          Tabular ouptus
```

### Model Evaluation and Animated Time Series
The BIOCLIM models are trained on presence data only, but we evaluated the fits using 300 random point locations across the extent of northern Europe defined in the study (100 sampled from each of the three periods). We carried out five-fold cross-validation for each model. For each run we determined a sensitivity threshold based by maximizing model sensitivity and specificity, and then derived the mean AUC, TSS and kappa for this threshold. The resulting evaluation statistcs can be found in the `mean_auc`, `mean_TSS` and `mean_kappa` colums of [this table](https://github.com/JesperBorrePedersen/Hamburgian_Culture_Climate_Analysis/blob/master/tables/model_evaluation.csv). For the final projections shown in the figures of the manuscript, we fitted three final models using the whole of each respective training dataset as outlined in the methods and derived a threshold based on sensitivity and specificity of each model. Again we calculated the AUC, TSS and kappa for these models, the corresponding values can be found in the likewise named colums of the [same table](https://github.com/JesperBorrePedersen/Hamburgian_Culture_Climate_Analysis/blob/master/tables/model_evaluation.csv). We then projected the thresholded suitability for each model across Northern Europe for the mean climate of the three time periods, as well as the climate at each 100 year time-step between 15.1 kyr BP and 13.0 kyr BP. The film strips for those animations are shown in Figure 5 of the manuscript. Animations of these time-series can be found [here](/figures/animations.md).

### Landscape Metrics
In addition to the visual analysis of the maps provided in the manuscript, we also calculated landscape metrics to describe the suitable climate-niche space for the three main model predictions. To achieve this, we first converted the raw suitability projections into a binary suitable / non-suitable landscape. For this we used the thershold values determined as outlined in the paragraph above. We then calculated the following landscape metrics using the `landscapemetrics` package in R (respective functions from the package in parentheses): proportion of total cells suitable (self-defined function, see code), number of patches (`lsm_l_np`), mean patch area (`lsm_c_area_mn`), standard deviation of the patch area (`lsm_c_area_sd`), mean distance to nearest neighbour (`lsm_p_enn`) and patch cohesion (`lsm_c_cohesion`). For simplicity we assumed that one cell of the CHELSA Trace21k raster was of a consistent (invariant) width and height of 1 km x 1 km. Please note that as the dataset comes in a geographic projection this is a very rough estimation and might lead to an overestimation of the real cell size. However, throughout the study extent the relative difference in width and height could assumed to be low, and should allow for a meaningful comparison among the areas and model projections shown.

### License
The content of this repository is licensed under a Creative Commons Attribution 4.0 International License: [CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)
