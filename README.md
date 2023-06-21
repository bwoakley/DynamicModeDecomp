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

## The algorithms

Suppose you have time series data $u_1 , u_2, \dots, u_N$ that evolves according to the linear map $u_{i+1} = A u_i$. If the data vectors $u_i$ are high dimensional, then so is the matrix $A$. 
If you are only given the time series data and the matrix $A$ is unknown, then approximating $A$ will prove difficult since it is high dimensional and not sparse. Instead of approximating the high dimensional $A$, the DMD algorithm computes a low dimensional projection $\widetilde{A}$ (see the write-up __Report.pdf__ for more details). The DMD algorithm in this form is due to Schmid 2010 (see [References](#references)).

The goal of this project is to extend the DMD algorithm to apply to time series data where the linearity hypothesis $u_{i+1} = A u_i$ fails. To this end, we propose several extensions of the DMD algorithm, they are: __learnB DMD__ , __learnAtilde DMD__, and __feedback DMD__. The write-up __Report.pdf__ is an introduction to these new algorithms and a measure of their improvement over the DMD algorithm by means of several detailed examples.

We now illustrate some of the methodology used in this research project.

## Data Cleaning

When learning the coefficients of the matrix $\widetilde{A}_i$ , we found that there was significant noise. To demonstrate that the entries of $\widetilde{A}_i$ evolved in a noisy way, we may consider the [Frobenius Norm](https://en.wikipedia.org/wiki/Matrix_norm) of differences: $$|| \widetilde{A}_{i+1}-\widetilde{A}_{i}||_F$$
We plot the evolution of this norm:


We can also demonstrate the noisy evolution of the entries of $\widetilde{A}_i$ by plotting these coefficients on top of each other:

As can be seen in the above plots, the evolution of the entries of $\widetilde{A}_i$ is quite noisy. Any model of this stochastic evolution would have a large variance. We may clean the data by altering one of the steps of the DMD algorithm to force the Singular Value Decomposition of $X_1$ to take its singular vectors from a codimension 1 submanifold (the right half of the 'plane'), so that the SVD is more unique. The relevant code, found in __DMD.m__, is:

	hi

This results in a more continuous evolutions of the coefficients of $\widetilde{A}_i$. We redraw the above plot after making this change. 
The evolution of the Frobenius Norm of the difference:

The evolution of the coefficients plotted on top of each other:

We recommend other researchers take this extra step before considering machine learning algorithms on $\widetilde{A}_i$ because the evolution is much smoother in the above two plots, making time series predictions more feasible. And so, we may pursue the proposed algorithm __learnAtilde DMD__. 

## References
- [ ] Schmid, P. (2010). Dynamic mode decomposition of numerical and experimental data. _Journal of Fluid Mechanics,_  _656_, 5-28. doi:[10.1017/S0022112010001217](https://doi.org/10.1017/S0022112010001217)
- [ ] Kutz, N. , Brunton, S. , Brunton, B. , and Proctor, J. (2016). Dynamic Mode Decomposition. _Society for Industrial and Applied Mathematics,_ Philadelphia, PA. doi:[10.1137/1.9781611974508](https://epubs.siam.org/doi/abs/10.1137/1.9781611974508)




## Project status
Complete.
Preliminary research results and conjectures to pursue in future work may be found in the document __DMDreport.pdf__ . 
