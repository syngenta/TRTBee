function save_results(path, fname, out)
% Short script to save TRT results
% Inputs:
% path - default path to save results
% fname - default file name
% out - cell array containing all relevant results from the TRT

% open ui box for user to select location:
[file, path] = uiputfile(strcat(path, fname));

if file ~=0
    % make sure there's a .csv still at the end
    [~, baseFileName, ext] = fileparts(file);
    if isempty(ext)
        ext = '.csv';
    end
    file = strcat(path, baseFileName, ext); % full path for file

    writecell(out, file);
end