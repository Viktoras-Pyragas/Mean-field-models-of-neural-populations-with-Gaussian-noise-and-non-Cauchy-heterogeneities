%----------------------------------------------------------------------------
% Normal approximation fixing the HWHM (Half width at half maximum)
% For normal distribution
%----------------------------------------------------------------------------
close all
clear 
d=1; % half width at half max (HWHM);
sig=d/(sqrt(2*log(2))); % variation (expressed through HWHM d);
my_color=0.6*[1 1 1];
mm=[1,2,4,6]; % changing the parameter of distribution;
N=500;
x=linspace(-5*d, 5*d,N);
nm=length(mm);
gr=zeros(nm,N);
for j=1:nm
    n=mm(j);
    [bet,C]=betn_Cn(n); % finding parameters of distribution for heterogeneity parameter n;
    % bet -> real positive coefficient;
    % C -> normalization constant;
    s=zeros(1,length(x));
    xx=ones(1,length(x));
    for k=1:n
        xx=bet*(x.^2).*xx/k;
        s=s+xx;
    end
    gm=C./(1+s); % distribution with heterogeneity parameter n; 
    gr(j,:)=gm; % collecting all distribution into one array;
end
gN=exp(-x.^2/(2*sig^2))/(sqrt(2*pi)*sig); % normal distribution;
figure
hold on;
%G=plot(x,gG,'-','Color',my_color,'LineWidth',3);
plot(x,gN,'-','Color',my_color,'LineWidth',3); % normal distribution;
plot(x,gr(1,:),'-.b','LineWidth',1) % distribution with heterogeneity parameter n=1;
plot(x,gr(2,:),'--r','LineWidth',1) % distribution with heterogeneity parameter n=2; 
plot(x,gr(3,:),':b','LineWidth',1)  % distribution with heterogeneity parameter n=4;
plot(x,gr(4,:),'-r','LineWidth',1)  % distribution with heterogeneity parameter n=6;

ylim([-0.01,0.5])
ax = gca;
ax.Layer='top';
lgd=legend({'normal','n=1','n=2','n=4','n=6'},'Location','northeast');
lgd.Position=[0.69    0.670    0.1643    0.2048];
              %lgd = legend('One','Two','Three','Four');
              lgd.FontSize = 8;
              legend boxoff;
xlabel('$\zeta$','Interpreter','latex','FontSize',10); % parameter of heterogeneity;
ylabel('$f_n(\zeta)$','Interpreter','latex','FontSize',10); % distribution function
set(gcf,'Units','centimeters')
set(gcf, 'PaperSize', [8.5,7]);
set(gcf,'Position',[2,2,8.5,7])
% set(gcf,'renderer','Painters');
%  print -depsc -tiff -r300 -vector normal_distr.eps



