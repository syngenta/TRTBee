function [par_out, par_best_fit, MLL, NRMSE, r2, SPPE] = setup_fit(substance, surv, concs, SEL, path, figSaveDir)
% Top-level script to run 



% function runs through fitting (or not fitting) and returns all output
clear global % remove previous global variables
global DATA glo % setup necessary global structure for BYOM
DATA{1} = surv; % survival data

DATA{2} = 0; % No other form of data

[par_out, par_best_fit, MLL, NRMSE, r2, SPPE] = fit_GUTS(concs, SEL, path); % run the BYOM script

if SEL == 1 % convert SEL into string for GUTS death mechanism
    deathMech = 'SD';
else
    deathMech = 'IT';
end
% create a copy of parameter space with a meaningful name
psName = strcat(figSaveDir, '\', substance, '_', deathMech, '_PS.mat'); 
matName = strcat('BYOM_GUTS_fit_2023_sel', num2str(SEL), '_PS.mat'); 