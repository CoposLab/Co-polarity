clear
close all
set(0,'DefaultFigureVisible','on')

% Load data
load('figure_data/figure3_data.mat')

% Define colors
green=[0,1,0.5];
grey=[0.4,0.6,0.7];
blue=[0.2,0.2,0.9];
gg=[linspace(green(1),grey(1),20);linspace(green(2),grey(2),20);linspace(green(3),grey(3),20)];
gb=[linspace(grey(1),blue(1),20);linspace(grey(2),blue(2),20);linspace(grey(3),blue(3),20)];
ggb=flip([gg gb]');

xvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};
yvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};

f1=figure(1);
f1.Position = [0 0 1100 700];
subplot(2,2,2)
h=heatmap(xvals,yvals,bindunbindSup);
h.Title='Probability Supracellular';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')

xvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};
yvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};


subplot(2,2,1)
h=heatmap(xvals,yvals,bindingSup);
h.Title='Probability Supracellular';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')


xvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};
yvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};

subplot(2,2,4)
h=heatmap(xvals,yvals,unbindingSup');
h.Title='Probability Supracellular';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')

yvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};
xvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};

f2=figure(2);
f2.Position = [0 0 1100 700];
subplot(2,2,3)
h=heatmap(xvals,yvals,bindunbindCoA');
h.Title='Probability Co-Aligned';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')

xvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};
yvals={'Rho binding times 1000','Rho binding times 10','No change','Rac binding times 10','Rac binding times 1000'};


subplot(2,2,1)
h=heatmap(xvals,yvals,bindingCoA');
h.Title='Probability Co-Aligned';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')


xvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};
yvals={'Rac unbinding times 1000','Rac unbinding times 10','No change','Rho unbinding times 10','Rho unbinding times 1000'};

subplot(2,2,4)
h=heatmap(xvals,yvals,unbindingCoA);
h.Title='Probability Co-Aligned';
h.YLabel='Cell 1';
h.XLabel='Cell 2';
colormap(ggb)
clim([0.2,1])
set(gcf,'color','w')