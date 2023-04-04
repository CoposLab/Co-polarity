function [Konx,Kony,Kfbx,Kfby,Koffx,Koffy] = spatialrates(kon,kfb,koff,a,b,s,beta,cond)
% Set spatially-dependent chemical reaction rate constants
%
steepness = 10;

assoc = 'on';
fb = 'on';
disassoc = 'on';
coupled = 'mixed';

switch coupled
    case {'on'}
        % association rate
        switch assoc
            case {'on'}
                Konx = 2*kon*(0.1+beta*a/max(a));
                Kony = 2*kon*(0.1+beta*b/max(b));
            case {'off'}
                Konx = kon*ones(length(s),1);
                Kony = kon*ones(length(s),1);
        end

        % feedback rate
        switch fb
            case {'on'}
                Kfbx = 2*kfb*(0.1+beta*a/max(a));
                Kfby = 2*kfb*(0.1+beta*b/max(b));
            case {'off'}
                Kfbx = kfb*ones(length(s),1);
                Kfby = kfb*ones(length(s),1);
        end

        % disassociation rate
        switch disassoc
            case {'on'}
                Koffx = 2*koff*(0.1+beta*(1-a/max(a)));
                Koffy = 2*koff*(0.1+beta*(1-b/max(b)));
                %Koffx = 2*koff*((1-tanh(steepness*(s-1)))/2 + 0.1)/1.1;
                %Koffy = 2*koff*((tanh(steepness*(s-1))+1)/2 + 0.1)/1.1;
            case {'off'}
                Koffx = koff*ones(length(s),1);
                Koffy = koff*ones(length(s),1);
        end
        
    case {'mixed'}
        % association rate
        switch assoc
            case {'on'}
                if (cond==0)
                    Konx = 2*kon*(0.1+beta*a/max(a));
                    Kony = 2*kon*(0.1+beta*b/max(b));
                else
                    Konx = 2*kon*((1-tanh(steepness*(s-5)))/2 + 0.1)/1.1;
                    Kony = 2*kon*((tanh(steepness*(s-5))+1)/2 + 0.1)/1.1;
                end
            case {'off'}
                Konx = kon*ones(length(s),1);
                Kony = kon*ones(length(s),1);
        end
        
        % feedback rate
        switch fb
            case {'on'}
                if (cond==0)
                    Kfbx = 2*kfb*(0.1+beta*a/max(a));
                    Kfby = 2*kfb*(0.1+beta*b/max(b));
                else
                    Kfbx = 2*kfb*((1-tanh(steepness*(s-5)))/2 + 0.1)/1.1;
                    Kfby = 2*kfb*((tanh(steepness*(s-5))+1)/2 + 0.1)/1.1;
                end
            case {'off'}
                Kfbx = kfb*ones(length(s),1);
                Kfby = kfb*ones(length(s),1);
        end
        
        % disassociation rate
        switch disassoc
            case {'on'}
                
                if (cond==0)
                    %Koffx = 2*koff*(0.1+beta*(2-a)/max(2-a));
                    %Koffy = 2*koff*(0.1+beta*(2-b)/max(2-b));
                    % attempted on 7/17 (otherwise Koffx, Koffy fall below 0)
                    Koffx = 2*koff*(0.1+beta*(1-a/max(a)));
                    Koffy = 2*koff*(0.1+beta*(1-b/max(b)));
                else
                    Koffy = 2*koff*((1-tanh(steepness*(s-5)))/2 + 0.1)/1.1;
                    Koffx = 2*koff*((tanh(steepness*(s-5))+1)/2 + 0.1)/1.1;
                end
            case {'off'}
                Koffx = koff*ones(length(s),1);
                Koffy = koff*ones(length(s),1);
        end
        
    case {'off'}
        % association rate
        switch assoc
            case {'on'}
                Konx = 2*kon*((tanh(steepness*(s-1))+1)/2 + 0.1)/1.1;
                Kony = 2*kon*((1-tanh(steepness*(s-1)))/2 + 0.1)/1.1;
            case {'off'}
                Konx = kon*ones(length(s),1);
                Kony = kon*ones(length(s),1);
        end

        % feedback rate
        switch fb
            case {'on'}
                Kfbx = 2*kfb*((tanh(steepness*(s-1))+1)/2 + 0.1)/1.1;
                Kfby = 2*kfb*((1-tanh(steepness*(s-1)))/2 + 0.1)/1.1;
            case {'off'}
                Kfbx = kfb*ones(length(s),1);
                Kfby = kfb*ones(length(s),1);
        end

        % disassociation rate
        switch disassoc
            case {'on'}
                Koffx = 2*koff*((1-tanh(steepness*(s-1)))/2 + 0.1)/1.1;
                Koffy = 2*koff*((tanh(steepness*(s-1))+1)/2 + 0.1)/1.1;
            case {'off'}
                Koffx = koff*ones(length(s),1);
                Koffy = koff*ones(length(s),1);
        end
end
end

