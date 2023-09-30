clear
close all
set(0,'DefaultFigureVisible','on')

load('racrho_binding_unbinding.mat')
bindunbind=vals1;

green=[0,1,0.5];
grey=[0.4,0.6,0.7];
blue=[0.2,0.2,0.9];
gg=[linspace(green(1),grey(1),20);linspace(green(2),grey(2),20);linspace(green(3),grey(3),20)];
gb=[linspace(grey(1),blue(1),20);linspace(grey(2),blue(2),20);linspace(grey(3),blue(3),20)];
ggb=flip([gg gb]');

xvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};
yvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};

figure(1)
h=heatmap(xvals,yvals,bindunbind);
h.Title='Proportion Leader/Follower';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')

load('racrho_binding.mat')

xvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};
yvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};

figure(2)
h=heatmap(xvals,yvals,bindingvals);
h.Title='Proportion Leader/Follower';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')

load('racrho_unbinding.mat')

xvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};
yvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};

figure(3)
h=heatmap(xvals,yvals,unbindingvals);
h.Title='Proportion Leader/Follower';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')