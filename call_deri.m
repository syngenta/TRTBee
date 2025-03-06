% %% BYOM function call_deri.m (calculates the model output)
% %
% %  Syntax: [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
% %
% % This function calls the explicit function(s) in <simplefun.html
% % simplefun.m> to calculate damage. Survival is subsequently calculated
% % from that in this function. As input, it gets:
% %
% % * _t_    the time vector
% % * _par_  the parameter structure
% % * _X0v_  a vector with initial states and one concentration (scenario number)
% % * _glo_ the structure with various types of information (used to be global)
% %
% % The output _Xout_ provides a matrix with time in rows, and states in
% % columns. This function calls <simplefun.html simplefun.m>. The optional
% % output _TE_ is not used in this package. 
% %
% % This function is for the reduced GUTS model. The external function
% % simplefun provides the scaled damage over time and the background
% % mortality. Mortality due to chemical stress is calculated in this
% % function as a form of 'output mapping'.
% %
% % * Author: Tjalling Jager
% % * Date: November 2021
% % * Web support: <http://www.debtox.info/byom.html>
% % * Back to index <walkthrough_guts.html>
% 
% %  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
% %  This source code is licensed under the MIT-style license found in the
% %  LICENSE.txt file in the root directory of BYOM. 
% 
% %% Start
% 
% function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
% 
% % These outputs need to be defined, even if they are not used
% Xout2    = []; % additional uni-variate output, not used in this example
% zvd      = []; % additional zero-variate output, not used in this example
% 
% %% Initial settings
% % This part organises a few things. Initial values can be determined by a
% % parameter (overwrite parts of _X0_), and zero-variate data can be
% % calculated. See the example BYOM files for more information. Note that
% % this version of call_deri only works with simplefun.
% 
% min_t = 500; % minimum length of time vector (affects ODE stepsize as well)
% 
% % Unpack the vector X0v, which is X0mat for one scenario
% X0 = X0v(2:end); % these are the intitial states for a scenario
% 
% %% Calculations
% % This part calls the explicit model in <simplefun.html simplefun.m>) to
% % calculate the output (the value of the state variables over time). There
% % is generally no need to modify this part. This version does NOT use an
% % ODE solver.
% 
% c       = X0v(1); % the concentration (or scenario number)
% t       = t(:);   % force t to be a row vector (needed when useode=0)
% t_rem   = t;      % remember the original time vector (as we will add to it)
% glo.timevar = 0;  % flag for time-varying exposure
% 
% if isfield(glo,'int_scen') && ismember(c,glo.int_scen) % is c in the scenario range global?
%     % then we have a time-varying concentration
%         
%     glo.timevar = 1; % tell simplefun that we have a time-varying treatment
%     % This is a time-saver as the isfield and ismember calls take some
%     % time. However, with simplefun the time gain is probably negligible.
%     
%     Tev = read_scen(-2,c,-1,glo); % use read_scen to derive actual exposure concentration
%     % the -2 lets read_scen know we need events, the -1 that this is also needed for splines!
%     min_t = max(min_t,length(Tev(Tev(:,1)<t(end),1))*2); 
%     % For very long exposure profiles (e.g., FOCUS profiles), we now
%     % automatically generate a larger time vector (twice the number of the
%     % relevant points in the scenario).
%     
%     % Add the time points from the exposure scenario to the time vector
%     % entered. This will usually add very little and cost a bit more
%     % calculation time. However, it is a safe idea to have those points in
%     % there. This is needed to get the exact same MLL as given by openGUTS.
%     % If this is slow, you can always comment it out.
%     t = unique([t;Tev(:,1)]);
%     
% end
% % For all GUTS cases, it is safe to use a long time vector. For the hazard
% % and proper model, we need to numerically integrate the hazard rate over
% % time. For IT as we need to find the maximum damage. For IT in reduced
% % models, we could make shortcuts for some scenarios, but since we don't
% % know how people may modify the models in derivatives or simplefun, it is
% % better to be general here. If the time span of the data set or simulation
% % is very long (and the number of points in Tev is not), consider
% % increasing the size of min_t.
% if length(t) < min_t % make sure there are at least min_t points
%     t = unique([t;(linspace(t(1),t(end),min_t))']);
% end
% 
% TE = +inf; % dummy for time of events
% 
% % to calculate damage, use an explicit function provided in simplefun
% Xout = simplefun(t,X0,par,c,glo);
% 
% 
% %% Output mapping
% % _Xout_ contains a row for each state variable. It can be mapped to the
% % data. If you need to transform the model values to match the data, do it
% % here. 
% %
% % This set of files is geared towards the use of the analytical solution
% % for the damage level in simplefun. Therefore, only the damage level and
% % background survival is available at this point, and the survival due to
% % the chemical is added for each death mechanism below.
% 
% mw   = par.mw(1);             % median of threshold distribution
% bw   = par.bw(1);             % killing rate
% Fs   = max(1+1e-6,par.Fs(1)); % fraction spread of the threshold distribution
% beta = log(39)/log(Fs);       % shape parameter for logistic from Fs
% S    = Xout(:,glo.locS);      % take background survival from the model
% Dw   = Xout(:,glo.locD);      % take the correct state variable for scaled damage
% 
% 
% switch glo.sel
%     case 1 % stochastic death
%         haz    = bw * max(0,Dw-mw);          % calculate hazard for each time point
%         cumhaz = cumtrapz(t,haz);            % integrate the hazard rate numerically
%         S      = S .* min(1,exp(-1*cumhaz)); % calculate survival probability, incl. background
%     
%     case 2 % individual tolerance, log-logistic distribution is used
%         mw = max(mw,1e-100); % make sure that the threshold is not exactly zero ...
%         % Make sure that Dw does not decrease over time (dead animals dont
%         % become alive). Uses built-in Matlab function cummax.
%         maxDw = cummax(Dw); % copy the vector Dw to maxDw and find cumulative maximum
%         S     = S .* (1 ./ (1+(maxDw/mw).^beta)); % survival probability
%         % the survival due to the chemical is multiplied with the background survival
% 
%     case 3 % mixed model
%         % This is a fast way for GUTS proper. Damage is calculated only
%         % once, as it is the same for all individuals. Survival for
%         % different NECs is calculated below. 
%         n          = 200; % number of slices from the threshold distribution
%         Fs2        = 999^(1/beta); % fraction spread for 99.9% of the distribution
%         z_range    = linspace(mw/(1.5*Fs2),mw*Fs2,n); % range of NECs to cover 99.9%
%         prob_range = ((beta/mw)*(z_range/mw).^(beta-1)) ./ ((1+(z_range/mw).^beta).^2); % pdf for the log-logistic (Wikipedia)
%         prob_range = prob_range / sum(prob_range); % normalise the densities to exactly one
%         S1         = zeros(length(t),1); % initialise the survival probability over time with zeros
%         for i = 1:n % run through the different thresholds
%             haz    = bw * max(0,Dw-z_range(i));  % calculate hazard for each NEC
%             cumhaz = cumtrapz(t,haz);            % integrate the hazard rate numerically
%             surv   = min(1,exp(-1*cumhaz));      % calculate survival probability
%             S1     = S1 + surv * prob_range(i);  % add to S1, weighted for this NECs prob density
%         end
%         S = S .* S1; % make sure to add background hazard
% end
% 
% 
% Xout(:,glo.locS) = S;           % replace correct state by newly calculated survival prob.
% [~,loct] = ismember(t_rem,t);   % find where the requested time points are in the long Xout
% Xout = Xout(loct,:);            % only keep the ones we asked for
% 
%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
%
% This function calls the explicit function(s) in <simplefun.html
% simplefun.m> to calculate damage. Survival is subsequently calculated
% from that in this function. As input, it gets:
%
% * _t_    the time vector
% * _par_  the parameter structure
% * _X0v_  a vector with initial states and one concentration (scenario number)
% * _glo_  the structure with various types of information (used to be global)
%
% The output _Xout_ provides a matrix with time in rows, and states in
% columns. This function calls <simplefun.html simplefun.m>. The optional
% outputs _TE_ (time of events caught by ODE solver) and _Xout_ (additional
% uni-variate data) are not used in this package. The _zvd_ for
% zero-variate data is also not used here. 
%
% This function is for the reduced GUTS model. The external function
% simplefun provides the scaled damage over time and the background
% mortality. Mortality due to chemical stress is calculated in this
% function as a form of 'output mapping'.
%
% * Author: Tjalling Jager
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, not used in this example
zvd      = []; % additional zero-variate output, not used in this example

%% Initial settings
% This part organises a few things. Initial values can be determined by a
% parameter (overwrite parts of _X0_), and zero-variate data can be
% calculated. See the example BYOM files for more information. Note that
% this version of call_deri only works with simplefun.

min_t = 500; % minimum length of time vector (affects ODE stepsize as well)

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario

%% Calculations
% This part calls the explicit model in <simplefun.html simplefun.m>) to
% calculate the output (the value of the state variables over time). There
% is generally no need to modify this part. This version does NOT use an
% ODE solver.

c       = X0v(1); % the concentration (or scenario number)
t       = t(:);   % force t to be a row vector (needed when useode=0)
t_rem   = t;      % remember the original time vector (as we will add to it)
glo.timevar = 0;  % flag for time-varying exposure

if isfield(glo,'int_scen') && ismember(c,glo.int_scen) % is c in the scenario range global?
    % then we have a time-varying concentration
        
    glo.timevar = 1; % tell simplefun that we have a time-varying treatment
    % This is a time-saver as the isfield and ismember calls take some
    % time. However, with simplefun the time gain is probably negligible.
    
    Tev = read_scen(-2,c,-1,glo); % use read_scen to derive actual exposure concentration
    % the -2 lets read_scen know we need events, the -1 that this is also needed for splines!
    min_t = max(min_t,length(Tev(Tev(:,1)<t(end),1))*2); 
    % For very long exposure profiles (e.g., FOCUS profiles), we now
    % automatically generate a larger time vector (twice the number of the
    % relevant points in the scenario).
    
    % Add the time points from the exposure scenario to the time vector
    % entered. This will usually add very little and cost a bit more
    % calculation time. However, it is a safe idea to have those points in
    % there. This is needed to get the exact same MLL as given by openGUTS.
    % If this is slow, you can always comment it out.
    t = unique([t;Tev(:,1)]);
    
end
% For all GUTS cases, it is safe to use a long time vector. For the hazard
% and proper model, we need to numerically integrate the hazard rate over
% time. For IT in reduced models, we could make shortcuts for some
% scenarios, but since we don't know how people may modify the models in
% derivatives or simplefun, it is better to be general here. If the time
% span of the data set or simulation is very long (and the number of points
% in Tev is not), consider increasing the size of min_t.
if length(t) < min_t % make sure there are at least min_t points
    t = unique([t;(linspace(t(1),t(end),min_t))']);
end

TE = +inf; % dummy for time of events

% to calculate damage, use an explicit function provided in simplefun
Xout = simplefun(t,X0,par,c,glo);

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 
%
% This set of files is geared towards the use of the analytical solution
% for the damage level in simplefun. Therefore, only the damage level and
% background survival is available at this point, and the survival due to
% the chemical is added for each death mechanism below.

mw   = par.mw(1);             % median of threshold distribution
bw   = par.bw(1);             % killing rate
Fs   = max(1+1e-6,par.Fs(1)); % fraction spread of the threshold distribution
beta = log(39)/log(Fs);       % shape parameter for logistic from Fs
S    = Xout(:,glo.locS);      % take background survival from the model
Dw   = Xout(:,glo.locD);      % take the correct state variable for scaled damage
        
switch glo.sel
    case 1 % stochastic death
        haz    = bw * max(0,Dw-mw);          % calculate hazard for each time point
        cumhaz = cumtrapz(t,haz);            % integrate the hazard rate numerically
        S      = S .* min(1,exp(-1*cumhaz)); % calculate survival probability, incl. background
    
    case 2 % individual tolerance, log-logistic distribution is used
        mw = max(mw,1e-100); % make sure that the threshold is not exactly zero ...
        % Make sure that Dw does not decrease over time (dead animals dont
        % become alive). Uses built-in Matlab function cummax.
        maxDw = cummax(Dw); % copy the vector Dw to maxDw and find cumulative maximum
        S     = S .* (1 ./ (1+(maxDw/mw).^beta)); % survival probability
        % the survival due to the chemical is multiplied with the background survival

    case 3 % mixed model
        % This is a fast way for GUTS proper. Damage is calculated only
        % once, as it is the same for all individuals. Survival for
        % different NECs is calculated below. 
        n          = 200; % number of slices from the threshold distribution
        Fs2        = 999^(1/beta); % fraction spread for 99.9% of the distribution
        z_range    = linspace(mw/(1.5*Fs2),mw*Fs2,n); % range of NECs to cover 99.9%
        prob_range = ((beta/mw)*(z_range/mw).^(beta-1)) ./ ((1+(z_range/mw).^beta).^2); % pdf for the log-logistic (Wikipedia)
        prob_range = prob_range / sum(prob_range); % normalise the densities to exactly one
        S1         = zeros(length(t),1); % initialise the survival probability over time with zeros
        for i = 1:n % run through the different thresholds
            haz    = bw * max(0,Dw-z_range(i));  % calculate hazard for each NEC
            cumhaz = cumtrapz(t,haz);            % integrate the hazard rate numerically
            surv   = min(1,exp(-1*cumhaz));      % calculate survival probability
            S1     = S1 + surv * prob_range(i);  % add to S1, weighted for this NECs prob density
        end
        S = S .* S1; % make sure to add background hazard
end

Xout(:,glo.locS) = S;           % replace correct state by newly calculated survival prob.
[~,loct] = ismember(t_rem,t);   % find where the requested time points are in the long Xout
Xout = Xout(loct,:);            % only keep the ones we asked for

