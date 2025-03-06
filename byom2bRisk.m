function output = byom2bRisk(surv, concs)
% Function to automatically convert and save survival data 
% into MOSAIC and Brisk format

% extract replicates from row names
trtmnts = surv(1,2:end);

concObs = size(concs, 1) - 1; % timepoints where concentration was measured, if 1 then assume constant concentration

survTimes = surv(2:end, 1); % survival observation times

replicate = [];
conc = [];
time = [];
Nsurv = [];

for i = 1:length(trtmnts)
    sVec = surv(2:end, i+1); % extract survival data
    indices = (isnan(sVec)==0); % get indices where this treatment was observed at those times
    
    relTimes = survTimes(indices); % times relevant to this replicate
    relSurv = sVec(indices);

    if concObs == 1
        % constant concentrations
        relConcs = concs(2, i+1)*ones(sum(indices), 1);
    else
        % in the most awkward cases, concentration and survival data might
        % be measured at different times, so we interpolate
        concTimes = concs(2:end, 1); %measured conc times
        interpConcs = interp1(concTimes,concs(2:end, i+1) , relTimes, 'linear', 'extrap');
        
        relConcs = interpConcs;

    end
    
    % Finally, build this into the table
    replicate = [replicate; i*ones(sum(indices), 1)];
    conc = [conc; relConcs];
    time = [time; relTimes];
    Nsurv = [Nsurv; relSurv];
end

% build final table
output = table(replicate, conc, time, Nsurv);