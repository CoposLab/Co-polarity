close all;
clear

% Loads each mat file, calculates times to polarize
% If save==1, saves results to two xls files
% First one ends with _Totals.xls: one row [avg_steps_c1,avg_steps_c2,avg_steps_total,avg_steps_samedir,avg_steps_lf,num_pol_c1,num_pol_c2,num_pol_yes,num_pol_lf]
% Second one ends with _AllSims: each row is one sim [c1 time, c2 time, yes time, LF time]

tic

save=1;
load_save_location='./simulation_results/timetopolarizeresults_signal/rhodownc1_rhoupc2/1000RhoOn_1000RhoOff';

num_runs=100;

all_polarize_times_c1=zeros(1,num_runs);
all_polarize_times_c2=zeros(1,num_runs);
all_polarize_times_yes=zeros(1,num_runs);
all_polarize_times_lf=zeros(1,num_runs);


for run=1:num_runs

    load(strcat(load_save_location, int2str(run), '.mat'));

    polarizedc1=0; %has cell1 polarized yet
    polc1_counter=0;
    polarizedc2=0;
    polc2_counter=0;
    pol_samedir=0; %have they polarized in the same direction yet
    samedir_counter=0; %how many timesteps in a row have they polarized in the same direction
    pol_lf=0;
    lf_counter=0;

    for t=1:Nt
        a1=a1all(:,t);
        a2=a2all(:,t);
        b1=b1all(:,t);
        b2=b2all(:,t);
        xC1=xC1all(:,t);
        yC1=yC1all(:,t);
        xC2=xC2all(:,t);
        yC2=yC2all(:,t);


        % Check if polarized
        a1New = a1;
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
        b1New = b1;
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
        a2New = a2;
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
        b2New = b2;
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


        if polarizedc1==0 && ~isempty(dirIndexa1) && ~isempty(dirIndexb1)
            [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
            medang1 = th(1,dirIndexa1);
            if polc1_counter==0 % this is the first time polarizing
                poldir_c1=medang1;
                polc1_counter=1;
            elseif abs(medang1-poldir_c1)<=pi/12 % not the first time polarizing, and same direction as it was before
                polc1_counter=polc1_counter+1;
            else % not the first time polarizing, but in a different direction now
                polc1_counter=1;
                poldir_c1=medang1;
            end
            if polc1_counter>=100
                polarizedc1=1;
                % sprintf('Cell 1 polarized after %d steps', t-100)
                all_polarize_times_c1(run)=t-100;
            end
        elseif polarizedc1==0 && polc1_counter>0 && (isempty(dirIndexa1) || isempty(dirIndexb1))
            polc1_counter=0;
        end
        if polarizedc2==0 && ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
            [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
            medang2 = th(1,dirIndexa2);
            if polc2_counter==0 % this is the first time polarizing
                poldir_c2=medang2;
                polc2_counter=1;
            elseif abs(medang2-poldir_c2)<=pi/12 % not the first time polarizing, and same direction as it was before
                polc2_counter=polc2_counter+1;
            else % not the first time polarizing, but in a different direction now
                polc2_counter=1;
                poldir_c2=medang2;
            end
            if polc2_counter>=100
                polarizedc2=1;
                % sprintf('Cell 2 polarized after %d steps', t-100)
                all_polarize_times_c2(run)=t-100;
            end
        elseif polarizedc2==0 && polc2_counter>0 && (isempty(dirIndexa2) || isempty(dirIndexb2))
            polc2_counter=0;
        end
        if ~isempty(dirIndexa1) && ~isempty(dirIndexb1) && ~isempty(dirIndexa2) && ~isempty(dirIndexb2) && pol_samedir==0
            [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
            medang1 = th(1,dirIndexa1);
            medang2 = th(1,dirIndexa2);
            angdiff = min(abs(medang1-medang2),abs(2*pi-abs(medang1-medang2)));
            angTolerance=pi/4;
            if angdiff < angTolerance
                samedir_counter=samedir_counter+1;
            else
                samedir_counter=0;
            end
            if samedir_counter>=100
                pol_samedir=1;
                % sprintf('Both polarized in same direction after %d steps', t-100)
                all_polarize_times_yes(run)=t-100;
            end
        else
            samedir_counter=0;
        end
        if ~isempty(dirIndexa1) && ~isempty(dirIndexb1) && ~isempty(dirIndexa2) && ~isempty(dirIndexb2) && pol_lf==0
            [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
            medang1 = th(1,dirIndexa1);
            medang2 = th(1,dirIndexa2);

            if xor(abs(medang1-3*pi/2)<pi/4, abs(medang2-pi/2)<pi/4)
                if lf_counter==0 && abs(medang1-3*pi/2)<pi/4 %c2 is the leader
                    lf_leader=2;
                    lf_counter=lf_counter+1;
                elseif lf_counter==0 && abs(medang2-pi/2)<pi/4 %c1 is the leader
                    lf_leader=1;
                    lf_counter=lf_counter+1;
                elseif lf_counter>0 && lf_leader==2 && abs(medang1-3*pi/2)<pi/4 %c2 is still the leader
                    lf_counter=lf_counter+1;
                elseif lf_counter>0 && lf_leader==1 && abs(medang2-pi/2)<pi/4 %c1 is still the leader
                    lf_counter=lf_counter+1;
                elseif lf_counter>0 && lf_leader==1 && abs(medang1-3*pi/2)<pi/4 %c2 is new leader
                    lf_counter=1;
                    lf_leader=2;
                elseif lf_counter>0 && lf_leader==2 && abs(medang2-pi/2)<pi/4 %c1 is new leader
                    lf_counter=1;
                    lf_leader=1;
                end

                if lf_counter>=100
                    pol_lf=1;
                    % sprintf('Polarized as LF after %d steps', t-100)
                    all_polarize_times_lf(run)=t-100;
                end
            else
                lf_counter=0;
            end
        else
            lf_counter=0;
        end


        if t==Nt-1
            % Calculate final difference in direction angles
            angTolerance=pi/4;
            strongAngTolerance=pi/5;
            if isempty(dirIndexa1) || isempty(dirIndexb1) %if c1 is not polarized at the end, set pol time to 0
                all_polarize_times_c1(run)=0;
            end
            if isempty(dirIndexa2) || isempty(dirIndexb2) %if c2 is not polarized at the end, set pol time to 0
                all_polarize_times_c2(run)=0;
            end

            if ~isempty(dirIndexa1) && ~isempty(dirIndexb1) && ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
                medang1 = th(1,dirIndexa1);
                medang2 = th(1,dirIndexa2);
                angdiff = min(abs(medang1-medang2),abs(2*pi-abs(medang1-medang2)));
                if ~(angdiff < angTolerance) % if they are not a "yes" at the end, set pol yes time to 0
                    all_polarize_times_yes(run)=0;
                end

                if ~xor(abs(medang1-3*pi/2)<pi/4, abs(medang2-pi/2)<pi/4) %if not LF at the end, set lf pol time to 0
                    all_polarize_times_lf(run)=0;
                end
            end
        end

    end
end

num_pol_c1=sum(all_polarize_times_c1~=0);
num_pol_c2=sum(all_polarize_times_c2~=0);
num_pol_yes=sum(all_polarize_times_yes~=0);
num_pol_lf=sum(all_polarize_times_lf~=0);

total_pol_time_c1=sum(all_polarize_times_c1);
total_pol_time_c2=sum(all_polarize_times_c2);
total_pol_time_yes=sum(all_polarize_times_yes);
total_pol_time_lf=sum(all_polarize_times_lf);

avg_steps_c1=0;
avg_steps_c2=0;
avg_steps_total=0;
avg_steps_samedir=0;
avg_steps_lf=0;
if num_pol_c1>0
    avg_steps_c1=total_pol_time_c1/num_pol_c1;
    sprintf('Average polarization time for cell 1 with %d runs: %d timesteps',num_pol_c1,avg_steps_c1)
end
if num_pol_c2>0
    avg_steps_c2=total_pol_time_c2/num_pol_c2;
    sprintf('Average polarization time for cell 2 with %d runs: %d timesteps',num_pol_c2,avg_steps_c2)
end
if num_pol_c1+num_pol_c2>0
    avg_steps_total=(total_pol_time_c1+total_pol_time_c2)/(num_pol_c1+num_pol_c2);
    sprintf('Average polarization time for one cell with %d runs: %d timesteps',num_pol_c1+num_pol_c2,avg_steps_total)
end
if num_pol_yes>0
    avg_steps_samedir=total_pol_time_yes/num_pol_yes;
    sprintf('Average yes polarization time with %d runs: %d timesteps',num_pol_yes,avg_steps_samedir)
else
    sprintf('No yesses')
end
if num_pol_lf>0
    avg_steps_lf=total_pol_time_lf/num_pol_lf;
    sprintf('Average LF polarization time with %d runs: %d timesteps',num_pol_lf,avg_steps_lf)
end

if save==1
    writematrix([avg_steps_c1,avg_steps_c2,avg_steps_total,avg_steps_samedir,avg_steps_lf,num_pol_c1,num_pol_c2,num_pol_yes,num_pol_lf],...
        strcat(load_save_location,'_Totals.xls'))
    writematrix([all_polarize_times_c1' all_polarize_times_c2' all_polarize_times_yes' all_polarize_times_lf'], ...
        strcat(load_save_location, '_AllSims.xls'))
end

toc