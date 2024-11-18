%------------------------------------------------------------------
% readme_files.txt
%------------------------------------------------------------------
% Here we work with theta neurons system subject to Gaussian noise;
% Here the heterogeneity is flat;
%------------------------------------------------------------------
% GN_dyn.m -> this is the main file; we slightly (step by step) change the system parameter
% (sig2 or del) and by using Newton iterations find the location of the fixed point and the 
% corresponding eigenvalues of the system jacobi matrix;
%------------------------------------------------------------------
% LM1_v.m -> here we find a fragment of the system Jacobi matrix at the
% given fixed point; Here we alse compute the corresponding fragment of
% the r.h.s. nonlinear function;
%------------------------------------------------------------------
% plot_GN_eigvs.m -> here we plot the results, the dependence of the real
% part of the leading eigenvalue vs a system parameter;
%------------------------------------------------------------------
% Lv1_v.m -> computing the vector-column L1(np); (corresponding to the np-th pole);
% These vectors are filling the last row of the Jacobi matrix; 
%------------------------------------------------------------------
% Bv1_v.m -> composing the row-vectors filling the lowest row of the
% Jacobi matrix;
%------------------------------------------------------------------
% FSt_v.m -> compute the r.h.s. function for the 
% relaxation equation: dS/dt=(-S+R)/taud;
%------------------------------------------------------------------
