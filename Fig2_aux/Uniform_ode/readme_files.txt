%--------------------------------------------------------------------------------
% readme_files.txt
%--------------------------------------------------------------------------------
% Here we integrate the ODEs of Fourier components with flat heterogeneity with 
% integer parameter n=2; (it can be changed);
%--------------------------------------------------------------------------------
% We expand the pdf into Fourier series and insert that expansion in the FPE,
% thus obtaining the infinity order ODE for Fourier exponents; This system is 
% truncated to a finite M-th order ODE (we take M=300); In these Eqs. the 
% Gaussian noise is involved; its amplitude is \sigma (in the code sig2=\sigma^2);
% For sig2=0 the ODE of Fourier components is attracted to the OA manifold;
%--------------------------------------------------------------------------------
% GN_dyn_vt.m -> we integrate the ODEs with n-th order flat heterogeneity;
% We integrate the truncated system of ODEs of Fourier components;
% Here we take into account the Gaussian noise (\sigma > 0);
%--------------------------------------------------------------------------------
% MFE_vt.m -> the r.h.s. of the ODEs of Fourier components for flat heterogeneity;
%--------------------------------------------------------------------------------
% plot_GN_v.m -> plotting the results of program GN_dyn_vt.m;
%--------------------------------------------------------------------------------
% GM_amt.mat -> data file wit initial conditions;
%--------------------------------------------------------------------------------
% GM_dyn.mat -> data file with dynamics of solution;
%--------------------------------------------------------------------------------
