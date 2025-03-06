function save_all_figs_app(dir, figName, closeFigs)
% Function saves all helpful figures from fitting the GUTS model
% Inputs:
% dir - directory in which the figures will be saved
% figName - name root for figure
% closeFigs - Boolean, 0 to leave figures open, 1 to close. Closing often
% enforced by the app


figHandles = findall(0,'Type','figure'); % Identify all current figures

% remove all instances which are not figures (e.g. the Matlab app window)
remove_list = [];
for i = 1:length(figHandles) 
    if isempty(figHandles(i).Number)
        % Then it is a waitbar or the app, don't save this
        remove_list = [remove_list, i];
    end
end
figHandles(remove_list) = []; % remove these


for i1 = 1:length(figHandles) 
    figNum   = num2str(get(figHandles(i1), 'Number'));
    % savefig(figHandles(i1), strcat(figName, figNum, '.fig'));
    fullPath = strcat(dir, "\", figName);
    saveas(figHandles(i1), strcat(fullPath, figNum, '.png'));
end
if closeFigs == 1
    close all; % close all figures
end