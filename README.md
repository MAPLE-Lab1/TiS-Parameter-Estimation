# TiS-Parameter-Estimation
A Tanks-in-Series Approach for Parameter Estimation of Lithium-ion Battery Models

# Overview
Parameter Estimation for an electrochemical model is generally challenging due to the nonlinear nature and computational complexity of the model equations.
In this work, we utilize the Tanks-in-Series (TiS) battery model for Li-ion batteries. TiS approach allows for faster parameter estimation with similar accuracy when compared to the original Pseudo two-dimentional (p2D) model.

More details of this work can be found in the article “A Tanks-in-Series Approach for Parameter Estimation of Lithium-ion Battery Models” published in Journal of Electrochemical Society (https://iopscience.iop.org/article/10.1149/1945-7111/ac6b5d).

# What's in the repository?
A demo code estimating two parameters using TiS model is provided in MATLAB 2020b platform. The TiS model simulation was performed using ODE15s syntax and genetic algorithm optimizer within the MATLAB was used to estimate parameters.

A similar demo was also created in python where TiS model simulation was performed using PyBaMM package (https://www.pybamm.org/).

