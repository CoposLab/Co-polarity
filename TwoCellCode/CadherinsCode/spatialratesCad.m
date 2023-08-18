% Set spatially-dependent kinetic rates
%
% Last updated: 7/10/2023
% Katie Levandosky

function [Konx,Kony,Konz,Kfbx,Kfby,Kfbz,Koffx,Koffy,Koffz] = spatialratesCad(kon,kfb,koff,a,b,s,beta,cond,bound,boundz)
steepness = 10;
amax = 10;
bmax = 10;

assoc = 'on';
fb = 'on';
disassoc = 'on';


% Re-formulate replace actin concentrations
aa = a.*(a<=amax) + amax.*(a>amax);
bb = b.*(b<=bmax) + bmax.*(b>bmax);

% association rate
switch assoc
    case {'on'}
        if (cond==0)
            % Konx = kon*(1+beta*aa);
            % Kony = kon*(1+beta*bb);
            Konx = kon*(1+beta(1)*aa);
            Konx(bound) = kon*(1+beta(2)*aa(bound));

            Kony = kon*(1+beta(1)*bb);
            Kony(bound) = kon*(1+beta(2)*bb(bound));

            % Konz = zeros(length(s),1);
            % Konz(boundz) = kon*ones(length(boundz),1);
            % Konz(boundz) = kon*(1+beta(2)*aa(boundz));
            Konz = kon*ones(length(s),1);
            % Konz = kon*(1+beta(2)*aa);
        else
            % periodic
            Konx = kon*(tanh(steepness*(s-1.875)) - tanh(steepness*(s-5.625)) + 0.2)/2.2;
            Kony = kon*(2 - tanh(steepness*(s-1.875)) + tanh(steepness*(s-5.625)) + 0.2)/2.2;

            % Konz = zeros(length(s),1);
            % Konz(boundz) = kon*ones(length(boundz),1);
            Konz = kon*ones(length(s),1);

            % no-flux
            %Konx = kon*((1-tanh(steepness*(s-5)))/2 + 0.1)/1.1;
            %Kony = kon*((tanh(steepness*(s-5))+1)/2 + 0.1)/1.1;
        end
    case {'off'}
        Konx = kon*ones(length(s),1);
        Kony = kon*ones(length(s),1);

        % Konz = zeros(length(s),1);
        % Konz(boundz) = kon*ones(length(boundz),1);
        Konz = kon*ones(length(s),1);
end

% feedback rate
switch fb
    case {'on'}
        if (cond==0)
            % Kfbx = kfb*(1+beta*aa);
            % Kfby = kfb*(1+beta*bb);
            Kfbx = kfb*(1+beta(1)*aa);
            Kfbx(bound) = kfb*(1+beta(2)*aa(bound));

            Kfby = kfb*(1+beta(1)*bb);
            Kfby(bound) = kfb*(1+beta(2)*bb(bound));

            Kfbz = zeros(length(s),1);
        else
            % periodic
            Kfbx = kfb*(tanh(steepness*(s-1.875)) - tanh(steepness*(s-5.625)) + 0.2)/2.2;
            Kfby = kfb*(2 - tanh(steepness*(s-1.875)) + tanh(steepness*(s-5.625)) + 0.2)/2.2;

            Kfbz = zeros(length(s),1);

            % no-flux
            %Kfbx = kfb*((1-tanh(steepness*(s-5)))/2 + 0.1)/1.1;
            %Kfby = kfb*((tanh(steepness*(s-5))+1)/2 + 0.1)/1.1;
        end
    case {'off'}
        Kfbx = kfb*ones(length(s),1);
        Kfby = kfb*ones(length(s),1);

        Kfbz = zeros(length(s),1);
end

% disassociation rate
switch disassoc
    case {'on'}
        if (cond==0)
            Koffx = koff*ones(length(s),1);
            Koffy = koff*ones(length(s),1);

            % Koffz = zeros(length(s),1);
            % Koffz(boundz) = koff*ones(length(boundz),1);
            Koffz = koff*ones(length(s),1);
        else
            % periodic
            Koffx = koff*(2 - tanh(steepness*(s-1.875)) + tanh(steepness*(s-5.625)) + 0.2)/2.2;
            Koffy = koff*(tanh(steepness*(s-1.875)) - tanh(steepness*(s-5.625)) + 0.2)/2.2;

            % Koffz = zeros(length(s),1);
            % Koffz(boundz) = koff*ones(length(boundz),1);
            Koffz = koff*ones(length(s),1);

            % no flux
            %Koffx = koff*((tanh(steepness*(s-5))+1)/2 + 0.1)/1.1;
            %Koffy = koff*((1-tanh(steepness*(s-5)))/2 + 0.1)/1.1;
        end
    case {'off'}
        Koffx = koff*ones(length(s),1);
        Koffy = koff*ones(length(s),1);

        % Koffz = zeros(length(s),1);
        % Koffz(boundz) = koff*ones(length(boundz),1);
        Koffz = koff*ones(length(s),1);
end
end

