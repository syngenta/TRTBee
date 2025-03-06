%%% set_options
% Script takes the option setting from the prelim_checks script in BYOM

%% Defines options structures
% First, the options structure for the various calculation routines are
% specified with default values. 

% Options for parameter optimisation (used in calc_optim)
opt_optim.fit      = 1; % fit the parameters (1), or don't (0)
opt_optim.it       = 1; % show iterations of the simplex optimisation (1, default) or not (0)
opt_optim.type     = 1; % optimisation method 1) simplex, 2) simulated annealing (experimental), 
                        % 3) swarm (experimental), 4) parameter-space explorer
opt_optim.simno    = 2; % for simplex: number of runs of fminsearch (starting again from previous best)
opt_optim.anntemp  = 1e-4; % for simulated annealing, stopping temperature
opt_optim.swno     = 200; % for swarm optimisation, number of particles
opt_optim.swit     = 20; % for swarm optimisation, number of iterations
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
opt_optim.ps_plots = 1; % when set to 1, makes intermediate plots to monitor progress of parameter-space explorer
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
opt_optim.ps_notitle = 0; % set to 1 to suppress plotting title on parameter-space plot
opt_optim.ps_slice = 0; % set to 1 to use slice sampler for main rounds (for use by Tjalling only!)
opt_optim.ps_dupl  = 1; % set to 1 to remove duplicates from sample in calc_parspace (slow for rough=0!)

% Options for plotting (used in calc_and_plot)
opt_plot.zvd     = 1; % turn on the plotting of zero-variate data (if defined)
opt_plot.sub     = 1; % switch for putting state variables as sub-plots into a figure (set to 1)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.cn      = 1; % if set to 1, connect model line to points (only for bw=1)
opt_plot.limax   = 0; % if set to 1, limit axes to the data set for each stage
opt_plot.sho     = 1; % set to 1 to show all scenarios, 0 to only plot model for scenarios with data
opt_plot.repls   = 1; % set to 1 to plot replicates, 0 to plot mean responses, 2 to make boxplots
opt_plot.notitle = 0; % set to 1 to suppress plotting titles on graphs with fits
opt_plot.outl    = 1; % set to 1 to identify outliers in the plot (points with weight zero; not for survival)
opt_plot.legsup  = 0; % set to 1 to suppress legends on fits
opt_plot.annot   = 0; % extra subplot in multiplot for fits: 1) box with parameter estimates, 2) overall legend
opt_plot.simulik = 1; % if set to 1, calculate the log-likelihood when fit=0 (in simulations)
opt_plot.statsup = []; % vector with states to suppress in plotting fits
opt_plot.transf  = 1; % set to 1 to calculate means and SEs including transformations (if repls=0)
opt_plot.y_zero  = 1; % set to 1 to force y-axis to start at zero
opt_plot.plt_rst = 0; % set to 1 to reset symbols/colours for each additional data set

% Options for local sensitivity analyses (used in local_sens)
opt_sens.type     = 1; % type of analysis 1) relative sensitivity 2) absolute sensitivity
opt_sens.state    = 0; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
opt_sens.step     = 0.05; % fraction change in each parameter's value
opt_sens.notitle  = 0; % set to 1 to suppress plotting titles on graphs
opt_sens.par_sel  = 0; % set to 1 to use fitted parameters only

% Options for asymptotical standard errors (used in calc_ase)
opt_ase.step      = 0.01; % relative stepsize

% Options for profiling likelihoods (used in calc_proflik)
opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 0; % number of sub-optimisations to perform to increase robustness
opt_prof.subrng   = 5; % maximum factor on parameters (both higher and lower) for sub-optimisations
opt_prof.brkprof  = 0; % when a better optimum is located, stop (1) or automatically refit (2)
opt_prof.verbose  = 1; % set to 0 to suppress output to screen, 2 to only suppress plots
opt_prof.subann   = 0; % set to 1 to use simulated annealing (followed by simplex) instead of suboptimisations
                       % setting this options means that subopt has no effect
opt_prof.saved    = 0; % set to 1 to plot/display saved profiles

% Options for the slice sampler (used in calc_slice)
opt_slice.thin     = 1; % thinning of the sample (keep one in every 'thin' samples)
opt_slice.burn     = 100; % number of burn-in samples (0 is no burn in)
opt_slice.slwidth  = 10; % initial width of the slice (Matlab default is 10)
opt_slice.alllog   = 0; % set to 1 to put all parameters on log-scale before taking the sample
opt_slice.testing  = 1; % make additional tests on the sample: moving average and autocorrelation
opt_slice.subplt   = 1; % create single plot with subplots (1) or many individual plots (0)

% Options for the calculation of confidence intervals (used in calc_conf and plot_guts)
opt_conf.type     = 3; % use values from slice sampler (1), likelihood region (2), or parspace explorer (3) to make intervals
opt_conf.samerr   = 0; % include sampling error in bounds for survival data (set samerr=1; requires statistics toolbox, Bayes only)
opt_conf.n_samerr = 20; % number of sub-sampling trials for each parameter set (if samerr=1)
opt_conf.sens     = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
opt_conf.state    = 0; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
opt_conf.lim_set  = 0; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_conf.n_lim    = 200; % size of limited set (likelihood-region and parspace only)
opt_conf.set_zero = {}; % parameter name (as cell array of strings) to set to zero for calc_conf (e.g., the background hazard)
opt_conf.use_par_out = 0; % set to 1 to use par as entered into the plotting function for CIs, rather than from saved set

% Options for the calculation of LCx/LPx values (used in calc_lcx_lim and calc_lpx_lim)
opt_lcx_lim.Feff  = 0.50; % effect level (>0 en <1), x/100 in LCx/LPx
opt_lcx_lim.plot  = 1; % set to 0 to NOT make a plot of LCx vs time or survival vs time at LPx
opt_lcx_lim.scen_plot = 1; % set to 0 to NOT make an extra plot of the exposure profile
opt_lcx_lim.scen_type = 4; % type of definition for the exposure scenario (as used in make_scen), if called with a file name
opt_lcx_lim.notitle   = 0; % set to 1 to suppress titles above ECx plots

% Options for the calculation of ECx and EPx values (used in calc_ecx, calc_epx and calc_epx_window)
opt_ecx.Feff      = [0.10 0.50]; % effect levels (>0 en <1), x/100 in ECx/EPx
opt_ecx.plot      = 1; % set to 0 to NOT make a plot of ECx vs time (or effects at different LPx)
opt_ecx.backhaz   = 'hb'; % parameter name (as string) to set to zero for LCx/LPx calculation to remove background mortality
opt_ecx.setzero   = {}; % parameter names (as cell array of strings) for extra parameters to be set to zero
opt_ecx.statsup   = []; % states to suppress from the calculations (e.g., locS)
opt_ecx.mf_range  = [1 10 100 1000]; % range for MFs to make plots with for calc_epx_window
opt_ecx.par_read  = 0; % when set to 1 read parameters from saved set, but do NOT make CIs
opt_ecx.batch_epx = 0; % when set to 1 use batch mode (no output to screen)
opt_ecx.batch_plt = 0; % when set to 1 create plots when using calc_epx_window_batch
opt_ecx.notitle   = 0; % set to 1 to suppress titles above ECx plots
opt_ecx.start_neg = 1; % set to 1 to start the moving window at minus window width
opt_ecx.Cminmax   = [1e-7 1e6]; % absolute boundaries for ECx (and its CI), [min>0,max]
opt_ecx.calc_int  = 0;  % integrate survival and repro into 1) RGR, or 2) survival-corrected repro (experimental!)
opt_ecx.rob_win   = 0; % set to 1 to use robust EPx calculation, rather than with fzero
opt_ecx.rob_rng   = [0.1:0.1:20 21:1:99 100:5:300]'; % range for calculation of robust EPx (smaller steps in interesting region)
opt_ecx.nomarker  = 0; % set to 1 to suppress markers for the ECx-time plot
opt_ecx.mf_crit   = 30; % MF trigger for flagging potentially critical profiles (EP10<mf_crit)
opt_ecx.prune_win = 0; % set to 1 to prune the windows to keep the interesting ones
opt_ecx.batch_eff = 0.1; % in batch mode, by default, check where effect exceeds 10%
opt_ecx.Tstep     = 1; % stepsize or resolution of the time window (default 1 day)
opt_ecx.saveall   = 0; % set to 1 to save all output (EPx for each element of the sample) in separate folder (only when making CIs with calc_epx_windows)
opt_ecx.plot_log  = 0; % set to 1 to plot ECx-time on log-log scale
opt_ecx.showprof  = 0;  % set to 1 to show exposure profile in plots at top row (calc_epx_window)
opt_ecx.lim_yax   = 1e4; % limit for y-axis with EPx values (calc_epx_window)
opt_ecx.id_sel    = [0 1 0]; % scenario to use from X0mat, scenario identifier, flag for ECx to use scenarios rather than concentrations
% Note on the id_sel option: the first value specifies the ID from which
% the initial values for the states need to be taken. The second value
% specifies the ID that the scenario needs to get, for calculating ECx or
% EPx. This determines which parameter set is used (when glo.names_sep is
% used). Third value is a flag to force calc_ecx to use a scenario for c to
% enter call_deri, rather than an actual concentration. This is needed when
% derivatives cannot handle c as a concentration, or when using
% glo.names_sep. The default setting takes the initial values from a
% control (ID=0), and calls it scenario 1, and uses c for actual
% concentrations in calc_ecx.

% Options for the joint likelihood confidence region (used in calc_likregion)
opt_likreg.skipprof = 0; % skip profiling step; use boundaries from saved likreg set (1) or profiling (2)
opt_likreg.chull    = 0; % set to 1 to plot convex hull that approximates 95% edges of the likelihood region
opt_likreg.axbnds   = 1; % bind axes on the bounds of the hyperbox (1), accepted sample (2), or inner region (3)
opt_likreg.burst    = 100; % number of random samples from parameter space taken every iteration
opt_likreg.lim_out  = 0; % set to 1 to sample from a smaller part of space (enough for forward predictions)
                         
% Options for plotting survival results for the GUTS packages (used in plot_guts) 
opt_guts.timeresp  = 1; % set 1 for multiplot with survival versus time for each treatment
opt_guts.doseresp  = 1; % set 1 for multiplot with survival versus concentration for each time point
opt_guts.timedeath = 1; % set 1 for multiplot with deaths-per-interval versus time for each treatment
opt_guts.single_dr = 0; % set 1 for single plot with survival versus concentration for each time point
opt_guts.single_td = 0; % set 1 for single plot with expected and observed deaths in each time interval
opt_guts.bw        = 0; % plot the single dose-response plot in black and white
opt_guts.cn        = 1; % use connecting lines between data and model curve when opt_guts.bw=1

% Options for the calculation of the intrinsic population growth rate
opt_pop.fscen   = [1 0.9 0.8]; % three scenarios with limited food
opt_pop.plt_fly = 1; % set to 1 for plotting on the fly

% Options for TKTD plotting
opt_tktd.repls    = 0; % plot individual replicates (1) or means (0)
opt_tktd.obspred  = 1; % plot predicted-observed plots (1) or not (0), (2) makes 1 plot for multiple data sets
opt_tktd.preds    = 0; % set to 1 to only plot predictions from X0mat without data
opt_tktd.addzero  = 0; % set to 1 to always add a concentration zero to X0mat
opt_tktd.max_exp  = 1; % set to 1 to maximise exposure/damage plots on exposure rather than damage
opt_tktd.notitle  = 0; % set to 1 to suppress titles above plots
opt_tktd.transf   = 1; % set to 1 to calculate means and SEs including transformations
opt_tktd.min      = 1; % set to 1 to show a dotted line for the control (lowest) treatment
opt_tktd.statsup  = []; % states to suppress from the plots (e.g., locS)
opt_tktd.flip     = 0; % set to 1 to flip row 1 and 2 of the plot around
opt_tktd.plotexp  = 1; % set to 1 to plot exposure profile as area in damage plots
opt_tktd.lim_data = 0; % set to 1 to limit axes to data
opt_tktd.set_ctrl = 0; % set to 1 to use separate control per data set for the dotted lines
opt_tktd.set_no   = 1; % set to 1 to specify 'Set i' in the titles of the subplots
opt_tktd.sppe     = 0; % set to 1 to calculate SPPEs (relative error at end of test)
opt_tktd.symsz    = 6; % size for the plotting symbol (default 6)
opt_tktd.pltbar   = []; % plot a bar at specified location, for selected scenarios
% Note: opt_tktd.pltbar should be a matrix, with scenario id's in the first
% row, start time for the bar in second, and stop time in third. This was
% used for the STARTRENS project to create a bar to indicate approximate
% hatching time.