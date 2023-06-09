# Temporal-Coherence-in-Electron-Microscopy
Programs used in the process of publication of "Temporal coherence envelope function of field emission in electron microscopy"

## Description
This project contains programmes that simulate and analyse images formed by Low Energy Electron Microscopy (LEEM). The principle and calculations are based on the following articles: 
* [Xuan Tan Nguyen, Michael S. Altman, Ultramicroscopy](https://doi.org/10.1016/j.ultramic.2023.113751)
* [K.M.Yu et al 2019 Ultramicroscopy 200, 160-168](https://doi.org/10.1016/j.ultramic.2019.01.015)
* [A B Pang et al 2009 J. Phys.: Condens. Matter 21 314006](https://doi.org/10.1088/0953-8984/21/31/314006)
* [S.M.Schramm et al 2012 Ultramicroscopy 115, 88-108](https://doi.org/10.1016/j.ultramic.2011.11.005)

A poster summarising this research can be found in the [link](https://drive.google.com/file/d/1OXvwEWZkI2jBYA_pVRFQ-NgTUow1s4UK/view?usp=sharing).

## Installation
The programmes are written in Python 3 with common libraries such as numpy, matplotlib. 

## LEEM Fourier Optics Simulation
### Resolution as a function of defocus
THe file "R(dz).py" contains the main programme to simulate the intensity profile of a given object function, at different defocus values specified by the user. Both standard (NAC) and aberration-correct (AC) LEEM parameters specified in [Nguyen, Altman](https://doi.org/10.1016/j.ultramic.2023.113751) are provided. The aperture angle and other parameters are also specified and can be altered by the user.

## Contributing
For contribution request, please email the author at xtnguyenaa@connect.ust.hk.

## Authors and acknowledgement
This is part of the UROP project on LEEM Fourier Optics at HKUST of the author under the supervision of professor [M.S. Altman](https://physics.ust.hk/eng/people_detail.php?pplcat=1&id=1). 

## Project status
Active
