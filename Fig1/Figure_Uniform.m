%----------------------------------------------------------------------------
% Uniform approximation fixing the HWHM (Half width at half maximum)
% For uniform distribution;
%----------------------------------------------------------------------------
close all
clear 
d=1; % half width at half max (HWHM);
my_color=0.6*[1 1 1]; % color of rectangular distribution;
mm=[1,2,6,10]; % changing the parameter of uniform distribution;
N=500;
x=linspace(-5*d, 5*d,N);
nm=length(mm);
gr=zeros(nm,N);
for j=1:nm
    n=mm(j); % integer parameter of uniform distribution;
Cr=ones(1,N)*n*sin(pi/2/n)/pi; % normalization parameter of the uniform distribution;
gr(j,:)=Cr./(1+x.^(2*n)); % uniform distribution with parameter m;
end
gR=zeros(1,N);
gR(1,200:300)=0.5; % rectangular distribution;
figure
hold on;
plot(x,gR,'-','Color',my_color,'LineWidth',3); % rectangular distribution;
plot(x,gr(1,:),'-.b','LineWidth',1) % uniform distribution with n=1
plot(x,gr(2,:),'--r','LineWidth',1) % uniform distribution with n=2
plot(x,gr(3,:),':b','LineWidth',1)  % uniform distribution with n=6
plot(x,gr(4,:),'-r','LineWidth',1)  % uniform distribution with n=10

ylim([-0.01,0.51])
ax = gca;
ax.Layer='top';
lgd=legend({'uniform','n=1','n=2','n=6','n=10'},'Location','northeast');
aa=lgd.Position
lgd.Position=[0.69    0.670    0.1643    0.2048];
              %lgd = legend('One','Two','Three','Four');
              lgd.FontSize = 8;
              legend boxoff;
xlabel('$\zeta$','Interpreter','latex','FontSize',10); % parameter of heterogeneity;
ylabel('$f_n(\zeta)$','Interpreter','latex','FontSize',10); % distribution function
set(gcf,'Units','centimeters')
set(gcf, 'PaperSize', [8.5,7]);
set(gcf,'Position',[2,2,8.5,7])
% print -depsc -tiff -r300 -vector uniform_distr.eps



