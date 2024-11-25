%-------------------------------------------------------------------------------
% readme_files.txt
%-------------------------------------------------------------------------------
% Here we compute and plot graphs, in which are compared post-transient dynamics
% of microscopic and macroscopic models; for macroscopic model we chose
% n=6 (Gaussian heterogeneity distribution), and for microscopic - the normal (n->infty) 
% heterogeneity distribution;
%-------------------------------------------------------------------------------
% This directory contais the following files:
%-------------------------------------------------------------------------------
% main_v.m -> main file; plots the whole graph;
%-------------------------------------------------------------------------------
% micro_Gauss_original_del_sig2.m -> integrates the microscopic equations for original parameters [del,sig2];
%-------------------------------------------------------------------------------
% micro_Gauss_changed_del.m -> integrates the microscopic equations for changed parameter del;
%-------------------------------------------------------------------------------
% micro_Gauss_changed_sig2.m -> integrates the microscopic equations for changed parameter sig2;
%-------------------------------------------------------------------------------
% mikro_original_del_sig2.mat -> results from microscopic equations for orginal parameters [del,sig2];
%-------------------------------------------------------------------------------
% mikro_changed_del.mat -> results from microscopic equations for changed parameter del;
%-------------------------------------------------------------------------------
% mikro_changed_sig2.mat -> results from microscopic equations for changed parameter sig2;
%-------------------------------------------------------------------------------
% mean_field_original_del_sig2.mat -> results from the mean field equations for orginal parameters [del,sig2];
%-------------------------------------------------------------------------------
% mean_field_changed_del.mat -> results from the mean field equations for changed parameter del;
%-------------------------------------------------------------------------------
% mean_field_changed_sig2.mat -> results from the mean field equations for changed parameter sig2;
%-------------------------------------------------------------------------------
% ksi_b.m -> computes the poles of heterogeneity ksi, and the expansion 
% coefficients bj;
%--------------------------------------------------------------------------------
% bet_n.m -> computes the coefficient beta for the number of poles n;
%--------------------------------------------------------------------------------
% betn_Cn.m -> computes the coefficients betan and Cn for a given 
% number of poles n;
%--------------------------------------------------------------------------------
% randpdf_my.m -> generates random numbers from a user defined distribution
%--------------------------------------------------------------------------------
% tight_subplot_v.m -> creates "subplot" axes with adjustable gaps and margins
%--------------------------------------------------------------------------------
