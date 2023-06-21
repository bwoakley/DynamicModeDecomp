# Dynamic Mode Decomposition

## Introduction
The following is a research project that I worked on as a Postdoctoral Scholar at the School of Mathematical and Statistical Sciences at Arizona State University. As a researcher in the field of fluid dynamics, I pursued a data driven point of view, using modern machine learning algorithms on the time series data of fluid flow states. In particular, my research focused on the Dynamic Mode Decomposition dimensionality reduction machine learning algorithm and its application to modeling complex time series data.

Fluid dynamics is an important field of study for many industries including: Aerospace, Energy (wind, hydroelectric, fossil fuels / oil), Automotive, Environmental, Pharmaceutical, Chemical and Petrochemical, Food and Beverage, etc. In general, this is implemented as a complex time series of data, where one must rely on mathematical theorems and algorithms to lead explainability and predictability of the phenomenon recognized in computational simulation models. 

## Description
The [Dynamic Mode Decomposition (DMD) algorithm](https://en.wikipedia.org/wiki/Dynamic_mode_decomposition) is a dimensionality reduction algorithm appearing in Schmid 2010 (see [References](#references)). Given simple dynamics (linear flows), it is known that one can use the DMD algorithm to estimate the future fluid velocity data. In this project, we extend this method to apply to weakly nonlinear flows. This project has applications to turbulent flow predictions such as weather forecasts. 






## Files
All code is written in Matlab (version R2022b: 9.13.0.2105380).

- [ ] DMDreport.pdf : A file detailing the project and results.
- [ ] MatlabCode\DMDflow.m : The main Matlab file which generates all images, videos, and results. Change various Boolean values from False to True to explore different sections of the code and produce some nice pictures.
- [ ]  MatlabCode\DMD.m : Helper function implementing DMD. Mostly borrowed from Kutz et. al. 2016 (see [References](#references)).
- [ ] MatlabCode\DMDflow.m :  Helper function reading the fluid flow data from the Cases folder.
- [ ] Cases : This folder contains the raw fluid flow data as binary files (bin01XXX). The data was generated in Fortran using a fourth order Runge–Kutta method. The Fortran file can be found at: [Diablo.f](https://www.damtp.cam.ac.uk/user/jrt51/files.html).

## Visuals



## References
- [ ] Schmid, P. (2010). Dynamic mode decomposition of numerical and experimental data. _Journal of Fluid Mechanics,_  _656_, 5-28. doi:[10.1017/S0022112010001217](https://doi.org/10.1017/S0022112010001217)
- [ ] Kutz, N. , Brunton, S. , Brunton, B. , and Proctor, J. (2016). Dynamic Mode Decomposition. _Society for Industrial and Applied Mathematics,_ Philadelphia, PA. doi:[10.1137/1.9781611974508](https://epubs.siam.org/doi/abs/10.1137/1.9781611974508)




## Project status
Complete.
Preliminary research results and conjectures to pursue in future work may be found in the document __DMDreport.pdf__ . 
