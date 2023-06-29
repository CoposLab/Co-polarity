
close all;

%Define colors
colorLength = 50;
white = [1,1,1];
red = [1,0,0];
blue = [0,0.5,1];
maroon = [0.4,0,0];
navy = [0,0,0.3];
yellow = [1,0.9,0];
darkyellow = [0.9,0.5,0];
whitered = [linspace(white(1),red(1),colorLength)',linspace(white(2),red(2),colorLength)',linspace(white(3),red(3),colorLength)'];
redmaroon = [linspace(red(1),maroon(1),colorLength)',linspace(red(2),maroon(2),colorLength)',linspace(red(3),maroon(3),colorLength)'];
whiteredmaroon = [whitered;redmaroon];
whiteblue = [linspace(white(1),blue(1),colorLength)',linspace(white(2),blue(2),colorLength)',linspace(white(3),blue(3),colorLength)'];
bluenavy = [linspace(blue(1),navy(1),colorLength)',linspace(blue(2),navy(2),colorLength)',linspace(blue(3),navy(3),colorLength)'];
whitebluenavy = [whiteblue; bluenavy];
myColors = [linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];
redblue = abs(whiteblue+whitered)./2;
redwhiteblue = [flip(whitered); whiteblue];
whiteyellow = [linspace(white(1),yellow(1),colorLength)',linspace(white(2),yellow(2),colorLength)',linspace(white(3),yellow(3),colorLength)'];
yellowdarkyellow = [linspace(yellow(1),darkyellow(1),colorLength)',linspace(yellow(2),darkyellow(2),colorLength)',linspace(yellow(3),darkyellow(3),colorLength)'];
whitedarkyellow = [whiteyellow;yellowdarkyellow];

% Define circles
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.93:0.01:1);
[Xcol,Ycol] = pol2cart(th,rad);
ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1]';
ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1]';
ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2]';
ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2]';
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
[Xsm,Ysm] = pol2cart(th,rad);
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.83:0.01:0.9);
[Xmid,Ymid] = pol2cart(th,rad);


% Concentric circles
% Cell 1
% figc1=figure(16);
figcells=figure(16);
surf(Xcol,Ycol,ZBranch1);
view(2)
colormap(whitebluenavy)
freezeColors;
freezeColors(colorbar);
clim([0,max(max(b1),max(b2))])
shading interp
hold on;
surf(Xmid,Ymid,ZBund1);
colormap(whitedarkyellow)
clim([0,max(max(a1),max(a2))])
freezeColors;
freezeColors(jicolorbar);
shading interp
grid off
% axis square
set(gca,'XTick',[], 'YTick', [])
% title('Cell 1 Combined: Blue=Branched, Red=Bundled')
% scatter(Xsm(boundC1),Ysm(boundC1),'black');
% hold off;
set(gca,'XColor','w')
set(gca,'YColor','w')
set(gcf,'color','w');



% Cell 2
surf(Xcol,Ycol-2,ZBranch2);
view(2)
colormap(whitebluenavy)
freezeColors;
freezeColors(colorbar);
clim([0,max(max(b1),max(b2))])
shading interp
surf(Xmid,Ymid-2,ZBund2);
colormap(whitedarkyellow)
freezeColors;
freezeColors(jicolorbar);
clim([0,max(max(a1),max(a2))])
shading interp
grid off
axis equal
set(gca,'XTick',[], 'YTick', [])
title('Blue=Branched, Yellow=Bundled')
% scatter(Xsm(boundC2),Ysm(boundC2)-2,'black');

allmax = max(max(max(a1),max(a2)),max(max(b1),max(b2)));

flipc2 = flip(boundC2);
for i=1:length(boundC1)
    plot3([Xsm(boundC1(i)) Xsm(flipc2(i))], [Ysm(boundC1(i)) Ysm(flipc2(i))-2],[allmax+1,allmax+1],'black')
end

hold off;
box off;
set(gca,'XColor','w')
set(gca,'YColor','w')
set(gcf,'color','w');



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
if ~isempty(dirIndex1)
    hold on;
    quiver(0,0,Xsm(dirIndex1),Ysm(dirIndex1),0,'color',[0 0 0])
    hold off;
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
if ~isempty(dirIndex2)
    hold on;
    quiver(0,-2,Xsm(dirIndex2),Ysm(dirIndex2),0,'color',[0 0 0])
    hold off;
end
