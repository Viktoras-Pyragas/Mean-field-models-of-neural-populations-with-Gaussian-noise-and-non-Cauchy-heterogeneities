Mean-field models of neural populations with Gaussian noise and non-Cauchy
heterogeneities

This is the code for the results and figures in our paper "Mean-field models of neural 
populations with Gaussian noise and non-Cauchy heterogeneities". 
They are written in Matlab, and require recent versions of MATLAB R2022b.

The directories:

/Fig1/ -> plotting uniform distribution;
/Fig4/ -> plotting normal distribution;

/Fig2/ and /Fig3/ -> plotting Hopf biffurcations and dynamics for uniform distribution;

/Fig5/ and /Fig6/  -> plotting Hopf biffurcations and dynamics for normal distribution;

The subdirectories:

/Fig2/ -> Hopf bifurcation for uniform distribution

/Fig3/ -> dynamics for uniform distribution;

/Fig2_aux/Uniform_Lin_Eqs/ -> Hopf bifurcation for uniform distribution;
Here we perform the Newton iterations with the linearized 
equations for Fourier components;

/Fig2_aux/Uniform_ode/ -> Here we integrate the nonlinear equations 
for Fourier components (for uniform distribution) in order to find the stable fixed point. 
This fixed point will be used in the Newton iterations;

The subdirectories:

/Fig5/ -> computing and plotting Hopf biffurcations for normal distribution;

/Fig6/ -> computing and plotting dynamics for normal distribution;

/Fig5aux/Gauss_Lin_Eqs/ -> Hopf bifurcation for normal distribution;
Here we perform the Newton iterations with the linearized 
equations for Fourier components;

/Fig5aux/Gauss_ode/ -> Here we integrate the nonlinear equations 
for Fourier components (for normal distribution) in order to find the stable fixed point. 
This fixed point will be used in the Newton iterations;
