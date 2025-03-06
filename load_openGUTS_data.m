function [byomData, vConcs] = load_openGUTS_data(fnm)
% Updated function to read in GUTS data according to the openGUTS format
% and output in BYOM format
% Inputs:
% fnm: File path for the selected GUTS data

%%

% Alternative method for older forms of matlab:
% [num] = readmatrix(fnm, "NumHeaderLines",1); % Skip over summary row

num = readtable(fnm); % Read in table
% first cell is A2, skip summary row
num = table2cell(num);
num = cell2mat(num);

% Find the gap between survival and treatment data
nanRow = find(isnan(num(1:end, 1)), 1) - 1; % nan rows interrupt survival and concentration data. Minus 1 to stop it appearing in surv

concsStart = find(~isnan(num((nanRow+1):end,1))); % nan ends when concentrations start
vConcs = num((nanRow+concsStart):end, :);

if size(vConcs,1) == 1
    % single row, constant concentration. Concentrations can be scenario
    % names
    surv = [[-1, num(end,2:end)]; num(1:nanRow,:)];
    vConcs = [[-1, num(end,2:end)]; vConcs];
else
    % variable exposure. This currently doesn't include BYOM options
    % for fitting background hazard on negative and solvent control
    numScen = size(vConcs, 2); % This is 1 larger than actual as it includes observation time
    scenarios = -1:(numScen-2); % -1 for BYOM format, -2 to correct counting "off by one"
    surv = [scenarios; num(1:nanRow,:)];
    vConcs = [scenarios; vConcs]; % Still has the minus 1, dont think this matters
end
byomData = surv; % output table