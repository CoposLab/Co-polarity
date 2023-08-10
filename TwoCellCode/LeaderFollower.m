set(0,'DefaultFigureVisible','off')

filename = './results/rhodown2/100thRhoOn';
maxnum=100;

angle=pi/4;

siglocation = 5*pi/4;
sigper = 0.4;


counteryes=0;
countersn=0;
counterno=0;
counter1np=0;
counter2np=0;
counterlf=0;
counterdist=0;

signal_polarized = 0;
countbund1=0;
countbund2=0;
tic
for i=51:maxnum
    close all
    fig = openfig(strcat(filename,'Scatter_',string(i),'.fig'));
    ohf = findobj(gcf);
    figaxes = findobj(ohf(1), 'Type', 'axes');

    cell1 = findobj(figaxes(2),'Type','line');
    xBundled1 = get(cell1(3),'XData');
    yBundled1 = get(cell1(3),'YData');
    xBranched1 = get(cell1(4),'XData');
    yBranched1 = get(cell1(4),'YData');

    cell2 = findobj(figaxes(1),'Type','line');
    xBundled2 = get(cell2(3),'XData');
    yBundled2 = get(cell2(3),'YData');
    xBranched2 = get(cell2(4),'XData');
    yBranched2 = get(cell2(4),'YData');

    % figure(2)
    % subplot(1,2,1)
    % hold on
    % plot(xBundled1,yBundled1)
    % plot(xBranched1,yBranched1)
    % hold off
    % subplot(1,2,2)
    % hold on
    % plot(xBundled2,yBundled2)
    % plot(xBranched2,yBranched2)
    % hold off

    a1=yBranched1;
    b1=yBundled1;
    a2=yBranched2;
    b2=yBundled2;

    [th,rad] = meshgrid((0:3.6:360)*pi/180,1);

    % Find median for cell 1
    a1New = yBranched1;
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
    a2New = yBranched2;
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

    % Calculate difference in direction angles
    angTolerance=pi/4;
    strongAngTolerance=pi/5;
    if isempty(dirIndex1) && isempty(dirIndex2)
        samedirection='2NP';
        angdiff=NaN;
        counter2np=counter2np+1;
    elseif isempty(dirIndex1) || isempty(dirIndex2)
        samedirection='1NP';
        angdiff=NaN;
        counter1np=counter1np+1;
    else
        medang1 = th(1,dirIndex1);
        medang2 = th(1,dirIndex2);
        angdiff = min(abs(medang1-medang2),abs(2*pi-abs(medang1-medang2)));
        if angdiff < angTolerance
            samedirection='yes';
            counteryes=counteryes+1;
        elseif (abs(medang1-3*pi/2)<strongAngTolerance && abs(medang2-pi/2)<strongAngTolerance)
            samedirection='strong no; collision';
            countersn=countersn+1;
        else
            samedirection='no';
            counterno=counterno+1;
        end

        if xor(abs(medang1-3*pi/2)<angle, abs(medang2-pi/2)<angle)
            counterlf=counterlf+1;
        end
    end

    if isempty(dirIndex1) && ~isempty(dirIndex2) && max(b1)>1
        medang2 = th(1,dirIndex2);
        if abs(medang2-pi/2)>angle
            counterdist=counterdist+1;
        end
    end
    if isempty(dirIndex2) && ~isempty(dirIndex1) && max(b2)>1
        medang1 = th(1,dirIndex1);
        if abs(medang1-3*pi/2)>angle
            counterdist=counterdist+1;
        end
    end
    if isempty(dirIndex1) && isempty(dirIndex2) && ((max(b1)>1 && max(a2)>1) || (max(b2)>1 && max(a1)>1))
        counterdist=counterdist+1;
    end

    if ~isempty(dirIndex2)
        medang2 = th(1,dirIndex2);
        if abs(medang2 - siglocation)<(sigper*2*pi)/2
            signal_polarized=signal_polarized+1;
        end
    end

    % if isempty(dirIndex1) && max(b1)>1 %bundled takes over cell 1
    %     countbund1=countbund1+1;
    % end
    % if isempty(dirIndex2) && max(b2)>1 %bundled takes over cell 2
    %     countbund2=countbund2+1;
    % end

end
toc
sprintf('%d yes, %d strong no, %d 1NP, %d 2NP\nNumber leader/follower: %d\nNumber dist. effort: %d\nPolarized: %d', counteryes, countersn, counter1np, counter2np, counterlf,counterdist,signal_polarized)

% sprintf('Number where bundled took over cell 1: %d\nNumber where bundled took over cell 2: %d', countbund1, countbund2)

% sprintf('Number leader/follower: %d',counterlf)
% sprintf('Number leader/follower np: %d', counterlf2)