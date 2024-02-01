% 2 cells

set(0,'DefaultFigureVisible','off')

filename = './simulation_results/results_celldifference/1_2Ka_1_2Kb_allC2/branchedbundledpromotion/0_9kb0_9kc2alpha50max';
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
both_sig_pol = 0;
countbund1=0;
countbund2=0;

countcw=0;
countccw=0;
tic
for i=1:maxnum
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
    if (a1New(1)~=0 && a1New(end)~=0)
        zeroInda1_1=find(a1New==0,1,'first');
        zeroInda2_1=find(a1New==0,1,'last');
        dirIndexa1=ceil((zeroInda1_1+zeroInda2_1)/2) - 50;
    else
        inda1_1=find(a1New~=0,1,'first');
        inda2_1=find(a1New~=0,1,'last');
        dirIndexa1=ceil((inda1_1+inda2_1)/2);
    end
    b1New = yBundled1;
    b1New(b1New<1)=0;
    if (b1New(1)~=0 && b1New(end)~=0)
        zeroIndb1_1=find(b1New==0,1,'first');
        zeroIndb2_1=find(b1New==0,1,'last');
        dirIndexb1=ceil((zeroIndb1_1+zeroIndb2_1)/2) - 50;
    else
        indb1_1=find(b1New~=0,1,'first');
        indb2_1=find(b1New~=0,1,'last');
        dirIndexb1=ceil((indb1_1+indb2_1)/2);
    end
    a2New = yBranched2;
    a2New(a2New<1)=0;
    if (a2New(1)~=0 && a2New(end)~=0)
        zeroInda1_2=find(a2New==0,1,'first');
        zeroInda2_2=find(a2New==0,1,'last');
        dirIndexa2=ceil((zeroInda1_2+zeroInda2_2)/2) - 50;
    else
        inda1_2=find(a2New~=0,1,'first');
        inda2_2=find(a2New~=0,1,'last');
        dirIndexa2=ceil((inda1_2+inda2_2)/2);
    end
    b2New = yBundled2;
    b2New(b2New<1)=0;
    if (b2New(1)~=0 && b2New(end)~=0)
        zeroIndb1_2=find(b2New==0,1,'first');
        zeroIndb2_2=find(b2New==0,1,'last');
        dirIndexb2=ceil((zeroIndb1_2+zeroIndb2_2)/2) - 50;
    else
        indb1_2=find(b2New~=0,1,'first');
        indb2_2=find(b2New~=0,1,'last');
        dirIndexb2=ceil((indb1_2+indb2_2)/2);
    end
    if dirIndexa1<1
        dirIndexa1=dirIndexa1+101;
    end
    if dirIndexb1<1
        dirIndexb1=dirIndexb1+101;
    end
    if dirIndexa2<1
        dirIndexa2=dirIndexa2+101;
    end
    if dirIndexb2<1
        dirIndexb2=dirIndexb2+101;
    end

    % Calculate difference in direction angles
    angTolerance=pi/4;
    strongAngTolerance=pi/5;
    if (isempty(dirIndexa1) || isempty(dirIndexb1)) && (isempty(dirIndexa2) || isempty(dirIndexb2))
        samedirection='2NP';
        angdiff=NaN;
        counter2np=counter2np+1;
    elseif (isempty(dirIndexa1) || isempty(dirIndexb1)) || (isempty(dirIndexa2) || isempty(dirIndexb2))
        samedirection='1NP';
        angdiff=NaN;
        counter1np=counter1np+1;
    else
        medang1 = th(1,dirIndexa1);
        medang2 = th(1,dirIndexa2);
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

        epsilon=0.01*2*pi;
        if xor(abs(medang1-3*pi/2)<angle, abs(medang2-pi/2)<angle) %leader/follower
            counterlf=counterlf+1;

            if abs(medang1-3*pi/2)<angle %if cell 2 is the leader
                if (medang2>=0 && medang2<(pi/2)-epsilon) || (medang2>(3*pi/2)+epsilon && medang2<=2*pi)
                    countccw=countccw+1;
                elseif medang2>(pi/2)+epsilon && medang2<(3*pi/2)-epsilon
                    countcw=countcw+1;
                end
            else %otherwise cell 1 is the leader
                if (medang1>=0 && medang1<(pi/2)-epsilon) || (medang1>(3*pi/2)+epsilon && medang1<=2*pi)
                    countcw=countcw+1;
                elseif medang1>(pi/2)+epsilon && medang1<(3*pi/2)-epsilon
                    countccw=countccw+1;
                end
            end
        end
    end

    if isempty(dirIndexa1) && ~isempty(dirIndexa2) && max(b1)>1
        medang2 = th(1,dirIndexa2);
        if abs(medang2-pi/2)>angle
            counterdist=counterdist+1;
        end
    end
    if isempty(dirIndexa2) && ~isempty(dirIndexa1) && max(b2)>1
        medang1 = th(1,dirIndexa1);
        if abs(medang1-3*pi/2)>angle
            counterdist=counterdist+1;
        end
    end
    if isempty(dirIndexa1) && isempty(dirIndexa2) && ((max(b1)>1 && max(a2)>1) || (max(b2)>1 && max(a1)>1))
        counterdist=counterdist+1;
    end

    if ~isempty(dirIndexa2)
        medang2 = th(1,dirIndexa2);
        if abs(medang2 - siglocation)<(sigper*2*pi)/2
            signal_polarized=signal_polarized+1;
        end
    end

    if ~isempty(dirIndexa1) && ~isempty(dirIndexa2)
        medang1 = th(1,dirIndexa1);
        medang2 = th(1,dirIndexa2);
        if abs(medang1 - siglocation)<(sigper*2*pi)/2 && abs(medang2 - siglocation)<(sigper*2*pi)/2
            both_sig_pol=both_sig_pol+1;
        end
    end

    % if isempty(dirIndexa1) && max(b1)>1 %bundled takes over cell 1
    %     countbund1=countbund1+1;
    % end
    % if isempty(dirIndexa2) && max(b2)>1 %bundled takes over cell 2
    %     countbund2=countbund2+1;
    % end

end
toc
sprintf(['%d yes, %d strong no, %d 1NP, %d 2NP\nNumber leader/follower: %d\n' ...
    'Number dist. effort: %d\nC2 polarized: %d\nBoth polarized: %d\n' ...
    'Number CW: %d, Number CCW: %d'], ...
    counteryes, countersn, counter1np, counter2np, counterlf,counterdist,signal_polarized,both_sig_pol,countcw,countccw)

% sprintf('Number where bundled took over cell 1: %d\nNumber where bundled took over cell 2: %d', countbund1, countbund2)

% sprintf('Number leader/follower: %d',counterlf)
% sprintf('Number leader/follower np: %d', counterlf2)