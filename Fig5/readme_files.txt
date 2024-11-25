%----------------------------------------------------------------------
% readme_files.txt
%----------------------------------------------------------------------
% We here take the case of normal heterogeneity distribution;
%----------------------------------------------------------------------
% Here we plot 4 bifurcation diagrams with supercritical Hopf
% bifurcations from two sources:
% 1) from averaged equations for truncated Fourier components; 
%    (with M=300 components);
% 2) from meanfield equations obtained from two-cumuliants
%    approach;
%----------------------------------------------------------------------
% There plotted 4 figures (a,b,c,d); Each figure corresponds to
% the following values of coupling strength J=(2,3,4,5);
% On each figure we plot Hopf curves for the following numbers of 
% heterogeneity: n=(1,2,4,6);
% On each figure the Hopf curves are plotted on the (\sigma^2,\Delta)
% plane;
%----------------------------------------------------------------------
% The gray solid lines are the Hopf curves computed from the mean-field
% equations derived from the two-cumulants approach;
% The colored lines are the Hopf curves computed from the truncated 
% equations for the Fourier components;
%----------------------------------------------------------------------
% The program contains the following functions:
%----------------------------------------------------------------------
% main_v.m -> the main function that organizes the plotting of the
% whole graph;
%----------------------------------------------------------------------
% interp_v2.m -> this function interpolates the results obtained from the
% truncated equations for the Fourier components;
% NOTE: this function first exchanges the arrays of del and sig2, and
% after interpolation the arrays are exchanged backwards again;
%----------------------------------------------------------------------
% myarrow_v.m -> plots an arrow; this function is used to plot
% an arrow from circle to the triangle, and from circle to the
% square;
%----------------------------------------------------------------------
