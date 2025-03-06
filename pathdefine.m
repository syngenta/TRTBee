%% BYOM function pathdefine.m (setting the path to the engine)
%
%  Syntax: pathdefine(varargin)
%
% A tiny piece of 'intelligent' code to add the required directories to the
% Matlab Path. Should work on UNIX, Mac and Windows, as long as the byom
% script from which this function is called is in a subdirectory of the
% BYOM directory. The addition to the path is temporary, and forgotten when
% Matlab is restarted.
% 
% In this version, the possibility is added to have TWO engine folders.
% This was done to allow use of the parallel computing toolbox of Matlab.
% By calling this function with a switch (in varargin), we can select the
% correct engine (and remove the other from the path). This requires that
% the additional engine is present in BYOM/engine_par.
% 
% Additionally, this version will remove any paths that contain
% openGUTS/engine. So no need to clean the path when switching from
% openGUTS to BYOM in the same Matlab session (openGUTS will be updated as
% well in the future to remove BYOM/engine path elements). 
%
% * Author: Tjalling Jager 
% * Date: April 2023
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>
% 
%  Copyright (c) 2012-2023, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function pathdefine(varargin)

if verLessThan('matlab','8.4')
    error('BYOM will not work for Matlab versions before 2014b!')
elseif verLessThan('matlab','9.1')
    warning('BYOM may produce errors for Matlab versions before 2016b!')
end

%% Search for BYOM folder
% This section checks whether BYOM is in the current directory path. If
% not, it gives an error. This should work on all OS types supported by
% Matlab, as it makes use of the Matlab built-in <filesep> for the file
% separator on the OS that is used.
   
curdir = pwd; % ask for path to current directory as string
k = strfind(curdir,'BYOM'); % find where BYOM is in the path to the current directory
if isempty(k) % if it is not in there, produce an error
    error('Make sure to place your scripts somewhere within the folder BYOM!')
end
warning('off','backtrace')
if length(k)>1 % then we have more folders with BYOM in the path!
    warning('Multiple folders with BYOM in the name have been found above this location; lowest is used!')
    k = k(end);
end
if ~strcmp(curdir(k-1:k(1)+4),[filesep,'BYOM',filesep])
    error('Cannot find the BYOM folder; make sure to not rename the BYOM folder!')
end

if isempty(varargin) % optional input argument, if it is missing ...
    usepar = 0; % set to 0 to NOT use separate engine folder for use with the parellel toolbox
else
    usepar = varargin{1}; % set to 1 to use separate engine folder for use with the parellel toolbox
end

% If the user requests to use parallel processing, that can only be done if
% the correct toolbox is installed, and when the correct BYOM version is
% present. If they are not, parallel processing is simply not used.
if usepar == 1
    v = ver; % take the version information
    if ~any(strcmp('Parallel Computing Toolbox',{v.Name}))
        % Note: the function <license> will only look at your *license* for
        % the toolbox, and *not* if it is actually installed!
        warning('You do not have the parallel computing toolbox installed, so the option for pathdefine is ignored!')
        usepar = 0;
    end
    if exist([curdir(1:k+3),filesep,'engine_par'],'dir') ~= 7 
        warning('Your BYOM platform does not contain the engine_par folder needed for parallel processing, so the option for pathdefine is ignored!')
        usepar = 0;
    end
end
warning('on','backtrace')

switch usepar % select the folders to use and the folders to remove
    case 0 % use standard engine
        add_dir = 'engine';
        rmv_dir = 'engine_par';
    case 1 % use engine_par folder
        add_dir = 'engine_par';
        rmv_dir = 'engine';
    otherwise
        error('Unknown option for pathdefine (use 0 or 1)')
end

if verLessThan('matlab','9.1') % this is to allow using old versions (before 2016b)!
    
    if ~isempty(strfind(path,['BYOM',filesep,add_dir,filesep])) %#ok<STREMP> correct engine is already in the path
        add_dir = []; % no need to add any directories
    end
    if isempty(strfind(path,['BYOM',filesep,rmv_dir,filesep])) %#ok<STREMP> wrong engine is already NOT in the path
        rmv_dir = []; % no need to remove any directories
    end
    
else % otherwise, use new recommended syntax
    
    if contains(path,['BYOM',filesep,add_dir,filesep]) % correct engine is already in the path
        add_dir = []; % no need to add any directories
    end
    if ~contains(path,['BYOM',filesep,rmv_dir,filesep]) % wrong engine is already NOT in the path
        rmv_dir = []; % no need to remove any directories
    end
    % Note: I look for the directory name with a <filesep> after it as,
    % otherwise, <contains> will also find 'engine' in 'engine_par'.
    %
    % Note: <contains> is the recommended syntax for Matlab. However,
    % <contains> was added in release 2016b, so this formulation will produce
    % an error in older versions. You can use strfind or findstr instead.
    
end

if isempty(add_dir) && isempty(rmv_dir)
    return % everything in the path is already set up correctly
end

%% Add engine folder to the path
% This section adds the correct engine folder, with its sub-folders, to the
% path. This should work on all OS types supported by Matlab, as it makes
% use of the Matlab built-in <filesep> for the file separator on the OS
% that is used.
   
% Calculate the path for BYOM, add and remove the correct engines to path
if ~isempty(rmv_dir)
    rmpath(genpath([curdir(1:k+3),filesep,rmv_dir])); % remove engine from the path, with its subdirectories
end
if ~isempty(add_dir)
    addpath(genpath([curdir(1:k+3),filesep,add_dir]),'-begin'); % add engine to the path, with its subdirectories
end

%% Remove any openGUTS elements from the path
% This will generally be superfluous, but there may be situations where you
% run openGUTS-Matlab and BYOM after each other. This section removes all
% elements with openGUTS/engine from the path.

if ~verLessThan('matlab','9.1') % skip this for really old versions!
    % older versions don't have <contains> but also not <split>
    
    if contains(path,['openGUTS',filesep,'engine']) % then user has previously ran the openGUTS Matlab version!
        path_cells = split(path,';'); % turn path into cell array
        for i = 1:length(path_cells)
            if contains(path_cells{i},['openGUTS',filesep,'engine'])
                rmpath(path_cells{i}); % remove openGUTS engine from the path
            end
        end
    end
    
end