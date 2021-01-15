# Hamburgian_Culture_Climate_Analysis

## Description

This repository contains the data and code for our paper: _Distribution modelling reveals the changing niches of pioneering Late Pleistocene populations in northern Europe._

## Content
1. [Authors](#Authors)
2. [Citation](#Citation)
3. [Acknowledgements](#Acknowledgements)
4. [Getting Started](#Getting-Started)
5. [Software](#Software)
6. [Folder Structure](#folder-structure)
7. [Animated Time Series](#animated-time-series)
8. [License](#License)

## Authors
Jesper Borre Pedersen, Jakob Johan Assmann, Signe Normand, Dirk Nikolaus Karger, Felix Riede.

Contact: jesper.borre@cas.au.dk

## Citation

Pedersen, Jesper B., Assmann, Jakob J., Normand, Signe, Karger, Dirk N. & Felix Riede (xxxx): _Distribution modelling reveals the changing niches of pioneering Late Pleistocene populations in northern Europe_.

## Acknowledgements
(from the manuscript)

Jesper B. Pedersen would like to thank the University of Aarhus for granting a PhD fellowship.

Felix Riede’s contribution to this paper is part of CLIOARCH, an ERC Consolidator Grant project that has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (grant agreement No. 817564).

## Getting Started
_Please note:_ This repository contains the code and the archaeological data required to replicate our analysis. However, the by-century temperature and preciptation rasters from the CHELSA Trace21k dataset are also needed to fully replicate the anlysis. These predictor data form the basis of the model fits and predictions. __The CHELSA Trace21k dataset is currenlty in review__ and not yet publicly available. We will update this section once the data has completed peer-review and has been released.  

The folder structure of this repostiroy and the software requirements for the anlysis are outlined below. 

All code required for the analysis and reproduction of the figures can be found in the `scripts/analysis.R` script. The code is structured using RStudio's section delination syntax to aid navigation within the script. 

### Software
The analysis was developed and carried out using R version 3.6.0 and the following packages: sf 0.9-3, raster 3.1-5, tidyverse 1.3.0, cowplot 1.0.0, rasterVis 0.47, colorspace 1.4-1, dismo 1.3-3, magick 2.5.2 and landscapemetrics 1.5.1 .

### Folder Structure

```
├── data            Archaeological data       
├── figures         Figure and animation outputs
├── scripts         The analysis script(s)
└── tables          Tabular ouptus
```

## Animated Time Series
The animated times series for the suitability predictions can be found [here](/figures/animations.md).

## License
The content of this repository is licensed under a Creative Commons Attribution 4.0 International License: [CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)
