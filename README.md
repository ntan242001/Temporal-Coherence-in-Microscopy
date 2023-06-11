# Temporal-Coherence-in-Electron-Microscopy
Programmes used in the process of publication of "[X. T. Nguyen, M. S. Altman, Temporal Coherence Envelope Function of Field Emission in Electron Microscopy, Ultramicroscopy 2023](https://doi.org/10.1016/j.ultramic.2023.113751))".

## Description
This project contains programmes that simulate and analyse images formed by Low Energy Electron Microscopy (LEEM). The principle and calculations are based on the article [Nguyen, Altman, Ultramicroscopy](https://doi.org/10.1016/j.ultramic.2023.113751), which considers the effect on image formation of partial temporal coherence from field emission sources used in most state-of-the-art electron microscopies.

Futher references include
* [K. M. Yu et. al. 2019 Ultramicroscopy 200, 160-168](https://doi.org/10.1016/j.ultramic.2019.01.015)
* [A. B. Pang et. al. 2009 J. Phys.: Condens. Matter 21 314006](https://doi.org/10.1088/0953-8984/21/31/314006)
* [S. M. Schramm et. al. 2012 Ultramicroscopy 115, 88-108](https://doi.org/10.1016/j.ultramic.2011.11.005)

A poster summarising this research can be found in the [link](https://drive.google.com/file/d/1OXvwEWZkI2jBYA_pVRFQ-NgTUow1s4UK/view?usp=sharing).

![image](https://github.com/ntan242001/Temporal-Coherence-in-Microscopy/assets/62791920/60108931-e0b0-404b-921b-54c8263c2b7a)

## Installation
The programmes are written in Python 3 with common libraries such as numpy, matplotlib. 

## LEEM Fourier Optics Simulation
### Resolution as a function of defocus
The file "R(dz).py" contains the main programme to simulate the intensity profile of a given object function, at different defocus values specified by the user. Both standard (NAC) and aberration-correct (AC) LEEM parameters specified in [Nguyen, Altman](https://doi.org/10.1016/j.ultramic.2023.113751) are provided. The aperture angle and other parameters are also specified and can be altered by the user.

### Figures
The codes used to plot the figures in the paper [Nguyen, Altman](https://doi.org/10.1016/j.ultramic.2023.113751) are available in the branch [Figures](https://github.com/ntan242001/Temporal-Coherence-in-Microscopy/tree/Figures).

## Contributing
For contribution request, please email the author at xtnguyenaa@connect.ust.hk.

## Acknowledgement
Support from the Hong Kong Research Grants Council (16304718) is gratefully acknowledged.

## Project status
Active
