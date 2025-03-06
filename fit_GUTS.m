function [par_out, par_best_fit, MLL, NRMSE, r2, SPPE] = fit_GUTS(concs, SEL, path)
% Run BYOM for inputted parameters:
% Inputs:
% concs - array of concentration values over time
% SEL - GUTS-RED death mechanism, either SD or IT
% fit_tox - fitting settings for BYOM
% path - path to save results

global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo

diary off           % turn of the diary function (if it is accidentaly on)

% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
% set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)

glo.basenm  = strcat(path, mfilename); % remember the filename for THIS file for the plots

glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

glo.fastrep = 0;


%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = unique(DATA{1}(1,2:end)); % scenarios (concentrations)
X0mat(2,:) = 1;                % initial survival probability
X0mat(3,:) = 0;                % initial scaled damage

% global parameters for GUTS purposes
glo.locS = 1; % location of survival probability in the state variable list
glo.locD = 2; % location of scaled damage in the state variable list

% start values and ranges are not used with parspace, but the fit flag and log setting are
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd = [1    1 1e-3  1e5 0];   % dominant rate constant, d-1
par.mw = [1    1    0  1e6 1];   % median threshold for survival (ug/L)
par.hb = [0.01 1    0    1 1];   % background hazard rate (1/d)
par.bw = [1    1 1e-6  1e6 0];   % killing rate (L/ug/d) (SD only)
par.Fs = [2    1    1  100 1];   % fraction spread of threshold distribution (IT only)
% NOTE: we use parspace for optimisation here, and an automated generation
% of search ranges. Therefore, the only aspect of the parameter structure
% _par_ that is used is the fit column.

if size(concs, 1) == 2
    X0mat(1,:) = concs(2, 2:end);
    DATA{1}(1,2:end) = concs(2, 2:end);
else
    Cw_type = 4; % 2 = block pulses, 4 = linear interp
    make_scen(Cw_type,concs);
end

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = ['scaled damage (',char(181),'g/L)'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'g/L']; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.
% opt_ecx.Feff = 0.5;
% opt_conf.type = 3; % use parspace explorer for CIs - calc_ecx


%% Calculations and plotting
% Here, the function is called that will do the calculation and the
% plotting. Options for the plotting can be set using opt_plot (see
% prelim_checks.m). Options for the optimisation routine can be set using
% opt_optim. The files in this directory always apply the analytical
% solution for damage in simplefun.


opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 1; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
% opt_plot.statsup = [2]; % vector with states to suppress in plotting fits
basenm_rem       = glo.basenm; % remember basename as we will modify it!
%%%%

%%%% Fitting controls
fit_tox = [-1 1]; % set by user (see table above)

opt_optim.type = 1; % optimisation method: 1) default simplex, 4) parspace explorer
par_out = automatic_runs_guts(fit_tox,par,[],SEL,[],opt_optim,opt_plot); % script to run the calculations and plot, automatically

par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% SEL     = [1 2]; % These are the death mechanisms that will be run automatically for fit_tox 1 and 2


fit_tox = [1 1]; % now an input to the function, meaning unchanged
opt_optim.type     = 4; % optimisation method 1) simplex, 4 parameter-space explorer

% Note: to use saved set for parameter-space explorer, use fit_tox(2)=0!
opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
opt_optim.ps_rough = 0; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS (2 for extra rough)
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
opt_optim.ps_plots = 0; % don't want any plots

conf_type = 3;

[par_out, MLL] = automatic_runs_guts(fit_tox,par,[],SEL,[],opt_optim,opt_plot); % script to run the calculations and plot, automatically
MLL = MLL(1);
opt_conf.type     = 3; % build CIs from parspace results
[h_coll, NRMSE, r2, SPPE] = plot_tktd_report(par_out,opt_tktd,opt_conf);
% Recalculate Goodness of fit so I can save it:
par_array = struct2cell(par_out);
par_array = cell2mat(par_array(1:(end-1))); 

if SEL == 1
    par_best_fit = par_array(1:(end-1),1)'; % Remove Fs from returned pars
else
    par_best_fit = par_array([1 2 3 5],1)'; % Remove bw from returned pars
end

par_array(par_array(:,5)==0,[1]) = log10(par_array(par_array(:,5)==0,[1])); % convert to log scale for relevant parameters

pfit = par_array(par_array(:,2)==1, 1); % only select parameters which were fitted

pmat = par_array;
[~,GoodFit,~] = transfer(pfit,pmat); %R2 matches printed output, but MLL doesn't. I favour the other method to extract MLL

