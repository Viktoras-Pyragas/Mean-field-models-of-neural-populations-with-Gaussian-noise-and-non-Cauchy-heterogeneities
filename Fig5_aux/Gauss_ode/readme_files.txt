%--------------------------------------------------------------------------------
% readme_ode_norm.txt
%--------------------------------------------------------------------------------
% Here we integrate the ODEs of Fourier components with normal heterogeneity with 
% integer parameter n=2; (it can be changed);
%--------------------------------------------------------------------------------
% We expand the pdf into Fourier series and insert that expansion in the FPE,
% thus obtaining the infinity order ODE for Fourier exponents; This system is 
% truncated to a finite M-th order ODE (we take M=300); In these Eqs. the 
% Gaussian noise is involved; its amplitude is \sigma (in the code sig2=\sigma^2);
% For sig2=0 the ODE of Fourier components is attracted to the OA manifold;
%--------------------------------------------------------------------------------
% GN_dyn_vt.m -> we integrate the ODEs with n-th order normal heterogeneity;
% We integrate the truncated system of ODEs of Fourier components;
% Here we take into account the Gaussian noise (\sigma > 0);
%--------------------------------------------------------------------------------
% plot_GN_v3.m -> breziame pr-mos GN_dyn_vt.m sk-mu rez-tus;
%--------------------------------------------------------------------------------
% MFE_vt.m -> the r.h.s. of the ODEs of Fourier components for flat heterogeneity;
%--------------------------------------------------------------------------------
% plot_GN_v.m -> plotting the results of program GN_dyn_vt.m;
%--------------------------------------------------------------------------------
% GM_amt.mat -> data file wit initial conditions;
%--------------------------------------------------------------------------------
% GM_dyn.mat -> data file with dynamics of solution;
%--------------------------------------------------------------------------------
% ksi_b.m -> computes the poles of heterogeneity ksi, and the expansion 
% coefficients bj;
%--------------------------------------------------------------------------------
% bet_n.m -> computes the coefficient beta for the number of poles n;
%--------------------------------------------------------------------------------
% betn_Cn.m -> computes the coefficients betan and Cn for a given 
% number of poles n;
%--------------------------------------------------------------------------------


