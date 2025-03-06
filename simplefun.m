% %% BYOM function simplefun.m (the model as explicit equations)
% %
% %  Syntax: Xout = simplefun(t,X0,par,c,glo)
% %
% % This function calculates the output of the reduced GUTS model system.
% % Note that the survival probability due to chemical stress is all
% % calculated in <call_deri.html call_deri.m>. This version can deal with
% % time-varying exposure concentrations, as long as the exposure profile is
% % specified using make_scen as type 2, 3 or 4 (smooth splining requires the
% % ODE version of the model in derivatives.m). If type 1 (interpolation) is
% % used for the exposure profile, make_scen will produce an error. I now
% % also added a make_scen type 4, which applies the analytical solution for
% % scaled damage under linear forcing sequentially (many thanks to Bob
% % Kooi). This allows to use an analytical solution for more complex
% % exposure scenarios as well.
% %
% % NOTE: watch out when the rate constant (_kd_ and _kc_) get close to each
% % other. In the analytical solution, there are divisions by the difference
% % between two rate constants, so they are not allowed to be equal. And,
% % there will be numerical problems when they get too close. I try to catch
% % these cases, and modify one or more of the rate constants a tiny bit.
% % However, it is good to be aware of this limitation (when in doubt: also
% % run the ODE version).
% %
% % As input, it gets:
% % 
% % * _t_   is the time vector
% % * _X0_  is a vector with the initial values for states
% % * _par_ is the parameter structure
% % * _c_   is the external concentration (or scenario number)
% % * _glo_ is the structure with information (normally global)
% %
% % Time _t_ is handed over as a vector, and scenario name _c_ as single
% % number, by <call_deri.html call_deri.m> (you do not have to use them in
% % this function). Output _Xout_ (as matrix) provides the output for each
% % state at each _t_.
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
% function Xout = simplefun(t,X0,par,c,glo)
% 
% %% Unpack initial states
% % The state variables enter this function in the vector _X0_. The initial
% % scaled damage level does not have to be zero. However, a non-zero value
% % would generally be meaningless (look at the full model for examples with
% % non-zero initial body residue). Furthermore, non-zero _Dw0_ will likely
% % lead to erroneous results for the calculation of LCx and LPx.
% 
% % S0  = X0(glo.locS); % survival probability at t=0
% Dw0 = X0(glo.locD); % scaled damage (referenced to external concentrations) at t=0
% 
% %% Unpack parameters
% % The parameters enter this function in the structure _par_. The names in
% % the structure are the same as those defined in the byom script file.
% % The 1 between parentheses is needed as each parameter has 5 associated
% % values.
% 
% kd   = par.kd(1);   % dominant rate constant
% % mw   = par.mw(1);   % median threshold (used in call_deri)
% % bw   = par.bw(1);   % killing rate (used in call_deri)
% % Fs   = par.Fs(1);   % fraction spread of NEC distribution (used in call_deri)
% hb   = par.hb(1);   % background hazard rate
% 
% %% Calculate the model output
% % This is the actual model, specified as explicit function(s):
% 
% S  = exp(-hb*t); % fill S with background mortality for survival
% % The calculation of the actual survival probability is done in call_deri
% % and only background survival is calculated here. 
% 
% % Do we need to do fast or slow kinetics?
% if isfield(glo,'fastslow') % does that field exist? than we may need to do something else
%     if glo.fastslow ~= 'o' % if it is not 'off' ...
%         c_v = c * ones(length(t),1); % if no exposure profile is specified, simply copy c over all time points
%         if glo.timevar == 1 % if we are warned that we have a time-varying concentration ...
%             c_v = read_scen(-3,c,t,glo); % use read_scen to derive actual exposure concentration vector
%             % the -3 lets read_scen know we are calling for fast/slow kinetics (and need a conc. vector)
%         end
%         if glo.fastslow == 'f'
%             Dw = c_v; % fast kinetics: damage equals external concentration
%             Dw(1) = Dw0; % make initial one Dw0 (especially needed for IT)
%         else
%             Dw = cumtrapz(t,c_v); % for slow kinetics, integrate over t
%         end
%         Xout(:,[glo.locS glo.locD]) = [S Dw]; % combine them into a matrix
%         return % return as we are done here
%     end
% end
% 
% % If regular kinetics ...
% Tev = [0 c]; % events setting: without anything else, assume it is constant 
% kc  = 0;     % assume no disappearance of the test chemical
% if glo.timevar == 1 % if we are warned that we have a time-varying concentration ...
%     [Tev,kc] = read_scen(-2,c,t,glo); % use read_scen to derive actual exposure concentration
%     % the -2 lets read_scen know we are calling from simplefun (and need events and kc)
% end
% 
% if sum(Tev(:,2)) == 0 && Dw0 == 0 % no need to calculate anything if there is only zero concentrations in Tev
%     Dw = zeros(length(t),1);
%     Xout(:,[glo.locS glo.locD]) = [S Dw]; % combine them into a matrix
%     return % return as we are done here
% end
% 
% % initialise internal concentrations and damage with NaNs
% Dw    = nan(length(t),1);
% Dw(1) = Dw0; % the first element is the starting concentration
% 
% diff_rel = 1e-5; % minimum relative difference between the rate constants
% 
% if size(Tev,2) == 2 % than we have a pulsed or static renewal scenario
% 
%     if abs(1-kd/kc) < diff_rel % the analytical solution does not allow the two rate constants
%         % to be exactly the same (or too close) ... when kc=0, the abs
%         % gives inf, so the part of the code below is not run.
%         kd = kd * (1+diff_rel); % increase kd tiny bit
%     end
% 
%     for i = 1:size(Tev,1) % run through all event periods
%         
%         a = kd * 1 * Tev(i,2); % modify the uptake flux to the new start concentration in this period
%         % compound parameter a is used in the solution below
%         
%         if i<size(Tev,1)
%             ind_t = (t>Tev(i,1) & t<=Tev(i+1,1)); % find the logical indices for the period
%         else
%             ind_t = (t>Tev(i,1)); % find the logical indices for the last period
%             % if the scenario is longer than the data, this is all zeros,
%             % which leads to an empty te, but this does not produce an
%             % error so it is fine.
%         end
%         te = t(ind_t) - Tev(i,1); % find the new part of the time vector, and make it start at zero again
%         
%         % analytical solution
%         Dw(ind_t) = (a/(kd-kc))*exp(-kc * te) + (Dw0 - (a/(kd-kc))) * exp(-kd * te);
%         
%         % find new starting values at exact moment of new period
%         if i < size(Tev,1) % don't do this for last period
%             te = Tev(i+1,1) - Tev(i,1); % start time for new period
%             Dw0 = (a/(kd-kc))*exp(-kc * te) + (Dw0 - (a/(kd-kc))) * exp(-kd * te);
%         end
%     end
%     
% else % we are using linear interpolation in a forcing series
%     
%     if t(end) > Tev(end,1) % do we ask for more points than in Tev?
%         Tev = [Tev; t(end) 0 0]; % add a dummy time point in Tev (from t)
%     end
% 
%     for i = 1:size(Tev,1)-1 % run through all events (not last one that we added)
%         
%         ind_t = (t>Tev(i,1) & t<=Tev(i+1,1)); % find the logical indices for the period
%         te = t(ind_t) - Tev(i,1); % find the new part of the time vector, and make it start at zero again
%         
%         % Thanks to Bob Kooi for using Maple to obtain the solution
%         a = Tev(i,2); % take as initial concentration the new one in this period
%         b = Tev(i,3); % take as slope the new one in this period
%         Dw(ind_t) = b * te + a - b / kd + exp(-kd * te) * (Dw0 - a + b / kd);
%         
%         % find new starting values at exact moment of new period
%         if i < size(Tev,1)-1 % don't do this for last period
%             te   = Tev(i+1,1) - Tev(i,1); % start time for new period
%             Dw0 = b * te + a - b / kd + exp(-kd * te) * (Dw0 - a + b / kd);
%         end
%     end
%     
% end
%    
% Xout(:,[glo.locS glo.locD]) = [S Dw]; % combine them into a matrix

%% BYOM function simplefun.m (the model as explicit equations)
%
%  Syntax: Xout = simplefun(t,X0,par,c,glo)
%
% This function calculates the output of the reduced GUTS model system.
% Note that the survival probability due to chemical stress is all
% calculated in <call_deri.html call_deri.m>. This version can deal with
% time-varying exposure concentrations, as long as the exposure profile is
% specified using make_scen as type 2, 3 or 4 (smooth splining requires the
% ODE version of the model in derivatives.m). If type 1 (interpolation) is
% used for the exposure profile, make_scen will produce an error. I now
% also added a make_scen type 4, which applies the analytical solution for
% scaled damage under linear forcing sequentially (many thanks to Bob
% Kooi). This allows to use an analytical solution for more complex
% exposure scenarios as well.
%
% NOTE: watch out when the rate constant (_kd_ and _kc_) get close to each
% other. In the analytical solution, there are divisions by the difference
% between two rate constants, so they are not allowed to be equal. And,
% there will be numerical problems when they get too close. I try to catch
% these cases, and modify one or more of the rate constants a tiny bit.
% However, it is good to be aware of this limitation (when in doubt: also
% run the ODE version).
%
% As input, it gets:
% 
% * _t_   is the time vector
% * _X0_  is a vector with the initial values for states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ is handed over as a vector, and scenario name _c_ as single
% number, by <call_deri.html call_deri.m> (you do not have to use them in
% this function). Output _Xout_ (as matrix) provides the output for each
% state at each _t_.
%
% * Author: Tjalling Jager 
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function Xout = simplefun(t,X0,par,c,glo)

%% Unpack initial states
% The state variables enter this function in the vector _X0_. The initial
% scaled damage level does not have to be zero. However, a non-zero value
% would generally be meaningless (look at the full model for examples with
% non-zero initial body residue). Furthermore, non-zero _Dw0_ will likely
% lead to erroneous results for the calculation of LCx and LPx.

% S0  = X0(glo.locS); % survival probability at t=0
Dw0 = X0(glo.locD); % scaled damage (referenced to external concentrations) at t=0

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

kd   = par.kd(1);   % dominant rate constant
% mw   = par.mw(1);   % median threshold (used in call_deri)
% bw   = par.bw(1);   % killing rate (used in call_deri)
% Fs   = par.Fs(1);   % fraction spread of NEC distribution (used in call_deri)
hb   = par.hb(1);   % background hazard rate

%% Calculate the model output
% This is the actual model, specified as explicit function(s):

S  = exp(-hb*t); % fill S with background mortality for survival
% The calculation of the actual survival probability is done in call_deri
% and only background survival is calculated here. 

% Do we need to do fast or slow kinetics?
if isfield(glo,'fastslow') % does that field exist? than we may need to do something else
    if glo.fastslow ~= 'o' % if it is not 'off' ...
        c_v = c * ones(length(t),1); % if no exposure profile is specified, simply copy c over all time points
        if glo.timevar == 1 % if we are warned that we have a time-varying concentration ...
            c_v = read_scen(-3,c,t,glo); % use read_scen to derive actual exposure concentration vector
            % the -3 lets read_scen know we are calling for fast/slow kinetics (and need a conc. vector)
        end
        if glo.fastslow == 'f'
            Dw = c_v; % fast kinetics: damage equals external concentration
            Dw(1) = Dw0; % make initial one Dw0 (especially needed for IT)
        else
            Dw = cumtrapz(t,c_v); % for slow kinetics, integrate over t
        end
        Xout(:,[glo.locS glo.locD]) = [S Dw]; % combine them into a matrix
        return % return as we are done here
    end
end

% If regular kinetics ...
Tev = [0 c]; % events setting: without anything else, assume it is constant 
kc  = 0;     % assume no disappearance of the test chemical
if glo.timevar == 1 % if we are warned that we have a time-varying concentration ...
    [Tev,kc] = read_scen(-2,c,t,glo); % use read_scen to derive actual exposure concentration
    % the -2 lets read_scen know we are calling from simplefun (and need events and kc)
end

if sum(Tev(:,2)) == 0 && Dw0 == 0 % no need to calculate anything if there is only zero concentrations in Tev
    Dw = zeros(length(t),1);
    Xout(:,[glo.locS glo.locD]) = [S Dw]; % combine them into a matrix
    return % return as we are done here
end

% initialise internal concentrations and damage with NaNs
Dw    = nan(length(t),1);
Dw(1) = Dw0; % the first element is the starting concentration

diff_rel = 1e-5; % minimum relative difference between the rate constants

if size(Tev,2) == 2 % than we have a pulsed or static renewal scenario

    if abs(1-kd/kc) < diff_rel % the analytical solution does not allow the two rate constants
        % to be exactly the same (or too close) ... when kc=0, the abs
        % gives inf, so the part of the code below is not run.
        kd = kd * (1+diff_rel); % increase kd tiny bit
    end

    for i = 1:size(Tev,1) % run through all event periods
        
        a = kd * 1 * Tev(i,2); % modify the uptake flux to the new start concentration in this period
        % compound parameter a is used in the solution below
        
        if i<size(Tev,1)
            ind_t = (t>Tev(i,1) & t<=Tev(i+1,1)); % find the logical indices for the period
        else
            ind_t = (t>Tev(i,1)); % find the logical indices for the last period
            % if the scenario is longer than the data, this is all zeros,
            % which leads to an empty te, but this does not produce an
            % error so it is fine.
        end
        te = t(ind_t) - Tev(i,1); % find the new part of the time vector, and make it start at zero again
        
        % analytical solution
        Dw(ind_t) = (a/(kd-kc))*exp(-kc * te) + (Dw0 - (a/(kd-kc))) * exp(-kd * te);
        
        % find new starting values at exact moment of new period
        if i < size(Tev,1) % don't do this for last period
            te = Tev(i+1,1) - Tev(i,1); % start time for new period
            Dw0 = (a/(kd-kc))*exp(-kc * te) + (Dw0 - (a/(kd-kc))) * exp(-kd * te);
        end
    end
    
else % we are using linear interpolation in a forcing series
    
    if t(end) > Tev(end,1) % do we ask for more points than in Tev?
        Tev = [Tev; t(end) 0 0]; % add a dummy time point in Tev (from t)
    end

    for i = 1:size(Tev,1)-1 % run through all events (not last one that we added)
        
        ind_t = (t>Tev(i,1) & t<=Tev(i+1,1)); % find the logical indices for the period
        te = t(ind_t) - Tev(i,1); % find the new part of the time vector, and make it start at zero again
        
        % Thanks to Bob Kooi for using Maple to obtain the solution
        a = Tev(i,2); % take as initial concentration the new one in this period
        b = Tev(i,3); % take as slope the new one in this period
        Dw(ind_t) = b * te + a - b / kd + exp(-kd * te) * (Dw0 - a + b / kd);
        
        % find new starting values at exact moment of new period
        if i < size(Tev,1)-1 % don't do this for last period
            te   = Tev(i+1,1) - Tev(i,1); % start time for new period
            Dw0 = b * te + a - b / kd + exp(-kd * te) * (Dw0 - a + b / kd);
        end
    end
    
end
   
Xout(:,[glo.locS glo.locD]) = [S Dw]; % combine them into a matrix