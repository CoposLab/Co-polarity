addpath('../../TwoCellCode/FigureAndMovieCode/freeze_colors')
set(0,'DefaultFigureVisible','on')
close all;
clear;
clc;

signal=0;
squished=0;
shifted=0;

setnum='1';
loadfile='./vid_matfiles/moving_cells_square/alt_racup_rhoup/1000RacOn_1000RhoOn';

load(strcat(loadfile,setnum));
t=2499;


allmax=max(max([a1all a2all a3all a4all b1all b2all b3all b4all]));

if shifted==0
    xshift1=0*xshift1;
    xshift2=0*xshift2;
    xshift3=0*xshift3;
    xshift4=0*xshift4;
    yshift1=0*yshift1;
    yshift2=0*yshift2;
    yshift3=0*yshift3;
    yshift4=0*yshift4;
end


a1=a1all(:,t);
a2=a2all(:,t);
b1=b1all(:,t);
b2=b2all(:,t);
a3=a3all(:,t);
a4=a4all(:,t);
b3=b3all(:,t);
b4=b4all(:,t);
L      = 10.0;
Na=101;



% Define circles
gapsize=0.01 * squished;
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.85:0.01:1);
[Xcol,Ycol] = pol2cart(th,rad);
Ycol1=Ycol;
Ycol2=Ycol;
Ycol3=Ycol;
Ycol4=Ycol;
Xcol1=Xcol;
Xcol2=Xcol;
Xcol3=Xcol;
Xcol4=Xcol;
if squished==1
    Ycol1(:,boundC1)=Ycol1(:,boundC1(1)*ones(1,length(boundC1)));
    Ycol2(:,boundC2)=Ycol2(:,boundC2(1)*ones(1,length(boundC2)));
end
Ycol2 = Ycol2 - 2*abs(max(max(Ycol2)))-gapsize;
Ycol3 = Ycol3 - 2*abs(max(max(Ycol3)))-gapsize;
Ycol4 = Ycol4 - 0*abs(max(max(Ycol4)))-gapsize;
Xcol3=Xcol3-2;
Xcol4=Xcol4-2;
ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1]';
ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1]';
ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2]';
ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2]';
ZBranch3 = [a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3]';
ZBund3 = [b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3]';
ZBranch4 = [a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4]';
ZBund4 = [b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4]';
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
[Xsm,Ysm] = pol2cart(th,rad);
Ysm1=Ysm;
Ysm2=Ysm;
if squished==1
    Ysm1(:,boundC1)=Ysm1(:,boundC1(1)*ones(1,length(boundC1)));
    Ysm2(:,boundC2)=Ysm2(:,boundC2(1)*ones(1,length(boundC2)));
end


sigper=0.40;
sigBound1 = (floor((Na-1)*3/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*3/8 + floor((Na-1)*sigper/2)))+1;
sigBound2 = (floor((Na-1)*5/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*5/8 + floor((Na-1)*sigper/2)))+1;

%Define colors
colorLength = 50;
white = [1,1,1];
darkyellow = [227/256,180/256,76/256];
yellow2 = [254/256,254/256,98/256];
pink = [211/256,95/256,183/256];
darkpink = [141/256,45/256,113/256];

whiteyellow2 = [linspace(white(1),yellow2(1),colorLength)',linspace(white(2),yellow2(2),colorLength)',linspace(white(3),yellow2(3),colorLength)'];
yellow2darkyellow = [linspace(yellow2(1),darkyellow(1),colorLength)',linspace(yellow2(2),darkyellow(2),colorLength)',linspace(yellow2(3),darkyellow(3),colorLength)'];
whitedarkyellow2 = [whiteyellow2;yellow2darkyellow];

whitepink = [linspace(white(1),pink(1),colorLength)',linspace(white(2),pink(2),colorLength)',linspace(white(3),pink(3),colorLength)'];
pinkdarkpink = [linspace(pink(1),darkpink(1),colorLength)',linspace(pink(2),darkpink(2),colorLength)',linspace(pink(3),darkpink(3),colorLength)'];
whitedarkpink = [whitepink;pinkdarkpink];


branchedColor = whitedarkpink;
bundledColor = whitedarkyellow2;
branchedColName = 'Pink';
bundledColName = 'Yellow';

% Concentric circles
% Cell 1
figcells=figure(1);
clf
hold on;
if squished==1
    plot3(cos(th(1,[1:boundC1(1),boundC1(end):end]))+xshift1(t+1),sin(th(1,[1:boundC1(1),boundC1(end):end]))+yshift1(t+1),ones(1,length(th(1,[1:boundC1(1),boundC1(end):end])))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
else
    plot3(cos(th(1,:))+xshift1(t+1),sin(th(1,:))+yshift1(t+1),ones(1,length(th))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
end
surf(Xcol1+xshift1(t+1),Ycol1+yshift1(t+1),ZBranch1,'AlphaData',ZBranch1,'FaceAlpha','interp','FaceColor','interp');
view(2)
colormap(branchedColor)
clim([0,allmax/2])
freezeColors;
shading interp
surf(Xcol1+xshift1(t+1),Ycol1+yshift1(t+1),ZBund1,'AlphaData',ZBund1,'FaceAlpha','interp','FaceColor','interp');
colormap(bundledColor)
freezeColors;
shading interp

% Cell 2
if squished==1
    plot3(cos(th(1,[1:boundC2(1),boundC2(end):end]))+xshift2(t+1),sin(th(1,[1:boundC2(1),boundC2(end):end]))-2*abs(max(max(Ycol2)))-gapsize+yshift2(t+1),ones(1,length(th(1,[1:boundC2(1),boundC2(end):end])))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
else
    plot3(cos(th(1,:))+xshift2(t+1),sin(th(1,:))-2+yshift2(t+1),ones(1,length(th))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
end
% plot3(cos(th(1,:))+xshift2(t+1),sin(th(1,:))-2+yshift2(t+1),ones(1,length(th))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
surf(Xcol2+xshift2(t+1),Ycol2+yshift2(t+1),ZBranch2,'AlphaData',ZBranch2,'FaceAlpha','interp','FaceColor','interp');
view(2)
colormap(branchedColor)
clim([0,allmax/2])
freezeColors;
shading interp
surf(Xcol2+xshift2(t+1),Ycol2+yshift2(t+1),ZBund2,'AlphaData',ZBund2,'FaceAlpha','interp','FaceColor','interp');
colormap(bundledColor)
clim([0,allmax/2])
freezeColors;
shading interp

% Cell 3
if squished==1
    plot3(cos(th(1,[1:boundC3(1),boundC3(end):end]))-2+xshift3(t+1),sin(th(1,[1:boundC3(1),boundC3(end):end]))-2*abs(max(max(Ycol3)))-gapsize+yshift3(t+1),ones(1,length(th(1,[1:boundC3(1),boundC3(end):end])))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
else
    plot3(cos(th(1,:))-2+xshift3(t+1),sin(th(1,:))-2+yshift3(t+1),ones(1,length(th))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
end
% plot3(cos(th(1,:))+xshift2(t+1),sin(th(1,:))-2+yshift2(t+1),ones(1,length(th))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
surf(Xcol3+xshift3(t+1),Ycol3+yshift3(t+1),ZBranch3,'AlphaData',ZBranch3,'FaceAlpha','interp','FaceColor','interp');
view(2)
colormap(branchedColor)
clim([0,allmax/2])
freezeColors;
shading interp
surf(Xcol3+xshift3(t+1),Ycol3+yshift3(t+1),ZBund3,'AlphaData',ZBund3,'FaceAlpha','interp','FaceColor','interp');
colormap(bundledColor)
clim([0,allmax/2])
freezeColors;
shading interp

% Cell 4
if squished==1
    plot3(cos(th(1,[1:boundC4(1),boundC4(end):end]))-2+xshift4(t+1),sin(th(1,[1:boundC4(1),boundC4(end):end]))-0*abs(max(max(Ycol4)))-gapsize+yshift4(t+1),ones(1,length(th(1,[1:boundC4(1),boundC4(end):end])))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
else
    plot3(cos(th(1,:))-2+xshift4(t+1),sin(th(1,:))+yshift4(t+1),ones(1,length(th))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
end
% plot3(cos(th(1,:))+xshift2(t+1),sin(th(1,:))-2+yshift2(t+1),ones(1,length(th))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
surf(Xcol4+xshift4(t+1),Ycol4+yshift4(t+1),ZBranch4,'AlphaData',ZBranch4,'FaceAlpha','interp','FaceColor','interp');
view(2)
colormap(branchedColor)
clim([0,allmax/2])
freezeColors;
cb=colorbar('Location','eastoutside');
freezeColors(cb);
cbpos=cb.Position;
set(cb,'Position',[cbpos(1)+2*cbpos(3),cbpos(2),cbpos(3),cbpos(4)/2])
set(cb,'TickLabels',[])
cbpos=cb.Position;
shading interp
surf(Xcol4+xshift4(t+1),Ycol4+yshift4(t+1),ZBund4,'AlphaData',ZBund4,'FaceAlpha','interp','FaceColor','interp');
colormap(bundledColor)
clim([0,allmax/2])
freezeColors;
jcb=jicolorbar;
freezeColors(jcb);
set(jcb,'Position',[cbpos(1)+cbpos(3),cbpos(2),cbpos(3),cbpos(4)])
shading interp
grid off
xlim([-5,2])
ylim([-5,3])
pbaspect([7 8 1])
% axis square

hold off;
box off;

bgColor=[1 1 1];
set(gca,'XTick',[],'YTick',[])
set(gca,'color',bgColor);
set(gcf,'color',bgColor);
set(gca,'XColor',bgColor)
set(gca,'YColor',bgColor)


% Find median for cell 1
a1New = a1;
a1New(a1New<1)=0;
if (a1New(1)~=0 && a1New(length(a1New))~=0)
    zeroInd1=find(a1New==0,1,'first');
    zeroInd2=find(a1New==0,1,'last');
    dirIndex1=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(a1New~=0,1,'first');
    ind2=find(a1New~=0,1,'last');
    dirIndex1=ceil((ind1+ind2)/2);
end
if dirIndex1<1
    dirIndex1=dirIndex1+101;
end

% Find median for cell 2
a2New = a2;
a2New(a2New<1)=0;
if (a2New(1)~=0 && a2New(length(a2New))~=0)
    zeroInd1=find(a2New==0,1,'first');
    zeroInd2=find(a2New==0,1,'last');
    dirIndex2=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(a2New~=0,1,'first');
    ind2=find(a2New~=0,1,'last');
    dirIndex2=ceil((ind1+ind2)/2);
end
if dirIndex2<1
    dirIndex2=dirIndex2+101;
end

% Find median for cell 3
a3New = a3;
a3New(a3New<1)=0;
if (a3New(1)~=0 && a3New(length(a3New))~=0)
    zeroInd1=find(a3New==0,1,'first');
    zeroInd2=find(a3New==0,1,'last');
    dirIndex3=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(a3New~=0,1,'first');
    ind2=find(a3New~=0,1,'last');
    dirIndex3=ceil((ind1+ind2)/2);
end
if dirIndex3<1
    dirIndex3=dirIndex3+101;
end

% Find median for cell 4
a4New = a4;
a4New(a4New<1)=0;
if (a4New(1)~=0 && a4New(length(a4New))~=0)
    zeroInd1=find(a4New==0,1,'first');
    zeroInd2=find(a4New==0,1,'last');
    dirIndex4=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(a4New~=0,1,'first');
    ind2=find(a4New~=0,1,'last');
    dirIndex4=ceil((ind1+ind2)/2);
end
if dirIndex4<1
    dirIndex4=dirIndex4+101;
end


% Plot arrows
if ~isempty(dirIndex1)
    hold on;
    quiver(0+xshift1(t+1),0+yshift1(t+1),Xsm(dirIndex1),Ysm(dirIndex1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7);
    hold off;
end
if ~isempty(dirIndex2)
    hold on;
    quiver(0+xshift2(t+1),-2+yshift2(t+1),Xsm(dirIndex2),Ysm(dirIndex2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7)
    hold off;
end
if ~isempty(dirIndex3)
    hold on;
    quiver(-2+xshift3(t+1),-2+yshift3(t+1),Xsm(dirIndex3),Ysm(dirIndex3),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7)
    hold off;
end
if ~isempty(dirIndex4)
    hold on;
    quiver(-2+xshift4(t+1),0+yshift4(t+1),Xsm(dirIndex4),Ysm(dirIndex4),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7)
    hold off;
end

% Plot signal
if signal==1
    [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
    [Xsig,Ysig] = pol2cart(th,rad);
    hold on;
    scatter(Xsig(sigBound2),Ysig(sigBound2)-2*abs(max(max(Ycol2)))-gapsize,'black','.')
    hold off;
end

% figure(2)
ohf = findobj(gcf);
figaxes = findobj(ohf(1), 'Type', 'axes');
set(figaxes(1),'Fontsize',15)
set(figaxes(2),'Fontsize',14)
camroll(90)