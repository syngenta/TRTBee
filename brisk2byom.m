function newName = brisk2byom(fnm)
% Function to automatically convert and save MOSAIC and Brisk format data 
% into BYOM format

% A = importdata('chlordan.txt');
A = importdata(fnm);
dt = A.data;

% extract replicates etc:
[repVals, id] = unique(dt(:,1)); % id gives indices in table to break up replicates
nReps = length(repVals); % number of replicates
allTimes = unique(dt(:,3)); % automatically sorts in ascending order

newDt = zeros(length(allTimes), nReps +1);
newDt(:,1) = allTimes;
newConc = zeros(length(allTimes), nReps +1);
newConc(:,1) = allTimes;

% for loop version. Might finesse later
for i = 1:nReps

    if i == nReps
        endId = length(dt(:,4));
    else
        endId = (id(i+1)-1);
    end
    % survivors:
    newDt(:,i+1) = dt(id(i):endId,4);
    % concentrations:
    newConc(:,i+1) = dt(id(i):endId,2);

end

% check for constant concentrations:
diffs = sum(diff(newConc(:,2:end))~=0); % look at differences in columns
if max(diffs) == 0
    % constant concs
    newConc = newConc(1,:);
end

% Combine all zeros into 1 treatment?

% stitch together the output with treatment names:
pre = 'T';
names = string();
for k = 1:nReps
    names = [names, strcat([pre,num2str(k,'%02d')])];
end
names(1) = []; % delete first entry which was used to set correct type

titles = ['Survival time [d]',  names]; 
titlesConc = ['Concentration time [d]',  names]; 


% arr = [titles; newDt; nan*ones(1,nReps+1)];
cLine = ['Concentration unit:', 'mug/bee/day', string(nan*ones(1, nReps-1))];
arr = [titles; newDt; nan*ones(1,nReps+1); cLine]; % vertcat
% 
arr = [arr; titlesConc; newConc];


% % add description
desc = "BLANK";
descLine = [desc, string(nan*ones(1, nReps))];
arr = [descLine;arr];
tbl = array2table(arr);

% rename and save
[~, baseFileName, ext] = fileparts(fnm);
newName = strcat(baseFileName,'_BYOM_format.txt');
writetable(tbl, newName, 'delimiter', '\t', 'WriteVariableNames', 0)

