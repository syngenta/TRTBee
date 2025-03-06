function run_fit_BYOM
% Top-level script to run within Matlab BYOM
% function runs through fitting (or not fitting) and returns all output
% clear global % remove previous global variables
% global DATA; % setup necessary global structure for BYOM



% Author     : Neil Sherborne
% Date       : June 2024

% Copyright (c) 2024-2034, Neil Sherborne, all rights reserved.
% This script, and the accompanying scripts in this package, is distributed
% under GNU GPL V3 License found in the LICENSE.txt file.

saveFigs = 1; % 1 to save figures, any other value won't save GUTS fit figures

fit_tox = [1, 1]; % Standard BYOM option for how to fit GUTS

[files, path] =uigetfile('*.txt','Select the INPUT DATA FILE(s)','MultiSelect','off');

if files ~=0 % zero returned if user cancels
    fnm2 = strcat(path, files);
    [~,substance]=fileparts(files); % remove file extension
    [surv, concs] = load_openGUTS_data(fnm2); % load data from the selected file
    DATA{1} = surv;
    meanDose = mean(concs(2:end,2:end), 1); % Calculate mean doses in case of time variable exposure
    if isempty(surv)
        error = strcat('No survival data extracted, check input file');
    end
end


figSaveDir = strcat(path, substance, "_GUTS_results"); % new directory to be created to save figures
% Check maximum mortality:
survEndProp = surv(end, 2:end)./surv(2, 2:end); % surviving proportion for all treatments
maxMort = (1 - min(survEndProp)); % maximum effect in the experiment
if maxMort <= 0.1
    mortFlag = 1;
else
    mortFlag = 0;
end
outFileName.Value = strcat(substance, "_TRT_analysis.csv");
maxDose = max(meanDose); % using nominal concentrations
%% Check for doses above 100 \mug/bee/d with less than 10% mortality
if (mortFlag && maxDose >= 100)
    Result = "Step 1: Low enough mortality at high dose (>100 \mug/bee) that TRT assessment is not needed";
    disp(Result);
    out = cell(1,2);
    out{1,1} = "Result: "; out{1,2} = Result;
    % add something to save results and return. Do this for the app as
    % well?
    % Can an LDD25 be calculated?
    
elseif maxMort >= 0.25
    %% Fit GUTS-RED-SD
    SEL = 1; % SD
    disp("Calibrating GUTS-RED-SD"); % Update progress in command window
    % Fit the GUTS model:
    [par_out_SD, par_best_fit_SD, MLL_SD, NRMSE_SD, r2_SD, SPPE_SD] = setup_fit(substance, surv, concs, SEL, path, figSaveDir);
    disp("Calculating SD LDDx values"); % Update progress window
    if saveFigs == 1
        % If we are saving figures...
        mkdir(figSaveDir) % Create folder for all figures and .mat files
        % Saving figure images...
        save_all_figs(figSaveDir, "SD_fit_result_figs", 1)
    end
    close all; % close all figures

    %% Calculate the LDDx values
    LDDxTimes = [10, 27, 182]; % Guidance needs LDD50 at 10, 27 and 182 days
    set_options
    opt_ecx.Feff = 0.50; % 50% effect for LDD50
    [LDD50_SD,LDD50_SD_lo,LDD50_SD_hi,Feff,ind_traits] = calc_ecx(par_out_SD,LDDxTimes,opt_ecx,opt_conf);
    % Output values into holders:
    LDD50_SD_10d = LDD50_SD{1}(LDDxTimes==10); LDD50_SD_10d_lo = LDD50_SD_lo{1}(LDDxTimes==10); LDD50_SD_10d_hi = LDD50_SD_hi{1}(LDDxTimes==10);
    LDD50_SD_27d = LDD50_SD{1}(LDDxTimes==27); LDD50_SD_27d_lo = LDD50_SD_lo{1}(LDDxTimes==27); LDD50_SD_27d_hi = LDD50_SD_hi{1}(LDDxTimes==27);
    LDD50_SD_182d = LDD50_SD{1}(LDDxTimes==182); LDD50_SD_182d_lo = LDD50_SD_lo{1}(LDDxTimes==182); LDD50_SD_182d_hi = LDD50_SD_hi{1}(LDDxTimes==182);

    %  LDD10s used for sublethal effects
    LDDxTimes = [27, 182];
    opt_ecx.Feff = 0.10;
    [LDD10_SD,LDD10_SD_lo,LDD10_SD_hi,Feff,ind_traits] = calc_ecx(par_out_SD,LDDxTimes,opt_ecx,opt_conf);
    LDD10_SD_27d = LDD10_SD{1}(1); LDD10_SD_27d_lo = LDD10_SD_lo{1}(1); LDD10_SD_27d_hi = LDD10_SD_hi{1}(1);
    LDD10_SD_182d = LDD10_SD{1}(2); LDD10_SD_182d_lo = LDD10_SD_lo{1}(2); LDD10_SD_182d_hi = LDD10_SD_hi{1}(2);

    % Working up to here :)
    %% Fit GUTS-RED-IT
    app.ProgressTextArea.Value = "Calibrating GUTS-RED-IT";
    SEL = 2; % IT
    [par_out_IT, par_best_fit_IT, MLL_IT, NRMSE_IT, r2_IT, SPPE_IT] = setup_fit(substance, surv, concs, SEL, path, figSaveDir);
    if saveFigs == 1
        % Saving figure images to the same directory
        save_all_figs(figSaveDir, "IT_fit_result_figs", 1)
    end

    %% Calculate the LDDx values
    app.ProgressTextArea.Value = "Calculating IT LDDx values";
    LDDxTimes = [10, 27, 182];
    set_options
    opt_ecx.Feff = 0.50;

    [LDD50_IT,LDD50_IT_lo,LDD50_IT_hi,Feff,ind_traits] = calc_ecx(par_out_IT,LDDxTimes,opt_ecx,opt_conf);
    LDD50_IT_10d = LDD50_IT{1}(LDDxTimes==10); LDD50_IT_10d_lo = LDD50_IT_lo{1}(LDDxTimes==10); LDD50_IT_10d_hi = LDD50_IT_hi{1}(LDDxTimes==10);
    LDD50_IT_27d = LDD50_IT{1}(LDDxTimes==27); LDD50_IT_27d_lo = LDD50_IT_lo{1}(LDDxTimes==27); LDD50_IT_27d_hi = LDD50_IT_hi{1}(LDDxTimes==27);
    LDD50_IT_182d = LDD50_IT{1}(LDDxTimes==182); LDD50_IT_182d_lo = LDD50_IT_lo{1}(LDDxTimes==182); LDD50_IT_182d_hi = LDD50_IT_hi{1}(LDDxTimes==182);

    LDDxTimes = [27, 182];
    opt_ecx.Feff = 0.10;
    [LDD10_IT,LDD10_IT_lo,LDD10_IT_hi,Feff,ind_traits] = calc_ecx(par_out_IT,LDDxTimes,opt_ecx,opt_conf);
    LDD10_IT_27d = LDD10_IT{1}(LDDxTimes==27); LDD10_IT_27d_lo = LDD10_IT_lo{1}(LDDxTimes==27); LDD10_IT_27d_hi = LDD10_IT_hi{1}(LDDxTimes==27);
    LDD10_IT_182d = LDD10_IT{1}(LDDxTimes==182); LDD10_IT_182d_lo = LDD10_IT_lo{1}(LDDxTimes==182); LDD10_IT_182d_hi = LDD10_IT_hi{1}(LDDxTimes==182);

    %% Check all quantitative errors against EFSA thresholds:
    minNRMSE = min(NRMSE_SD, NRMSE_IT);
    maxSPPE = max(max(abs(SPPE_SD{1})), max(abs(SPPE_IT{1})));

    if (minNRMSE > 0.5) || (maxSPPE > 50)
        % Both GUTS models failed to fit the data
        % Only option is to fit DRC
        % Worst case slope of -2 applies
        [LDDs, slope, LD50_emp] = haber_slope_fit(surv, concs, 0.9); % 0.9 for 10% effect

        out = cell(5); % setup output cell array
        Result = "Step 4: Neither GUTS-RED-SD, GUTS-RED-IT can sufficiently fit the data";
        Result2 = "Default slope of -2 assumed.";
        out{1,1} = "Result: "; out{1,2} = Result; out{2,2} = Result2;

        out{4,1} = "DRC slope"; out{4,2} = "DRC 10d LD50 (inflection point)";
        out{4,3} = "27d LDD10"; out{4,4} = "182d LDD10";
        out{5,1} = slope; out{5,2} = LD50_emp;
        out{5,3} = LDDs(1); out{5,4} = LDDs(2); % 27 and 182d LDD10

    else
        % At least one model fits sufficiently well
        if NRMSE_SD <= NRMSE_IT
            preferred_model = "SD";
        else
            preferred_model = "IT";
        end

        % Assess TRT according to Eq. (22) of bee guidance
        SD_TRT_assess = (LDD50_SD_27d < LDD50_SD_10d/2.7); % Boolean, true if TRT
        IT_TRT_assess = (LDD50_IT_27d < LDD50_IT_10d/2.7); % 1 if TRT



        %% print the info that will be in the output anyway
        % SD:
        out = cell(7,8);
        out{1,3} = "best fit:"; out{1,4} = preferred_model;
        out{3,1} = "SD parameters (kd, mw, bw, hb)"; out{3,2} = "SD NRMSE"; out{3,3} = "10d LDD50"; out{3,4} = "27d LDD50"; out{3,5} = "182d LDD50"; out{3,6} =  "10d LDD50/2.7"; out{3,7} = "27d LDD10"; out{3,8} = "182d LDD10";
        out{4,1} = par_best_fit_SD([1,2,4,3]); % hb last
        out{4,2} = NRMSE_SD; out{4,3} = strcat(num2str(LDD50_SD_10d), " (", num2str(LDD50_SD_10d_lo), ", ", num2str(LDD50_SD_10d_hi), ")");
        out{4,4} = [strcat(num2str(LDD50_SD_27d), " (", num2str(LDD50_SD_27d_lo), ", ", num2str(LDD50_SD_27d_hi), ")")];
        out{4,5} = [strcat(num2str(LDD50_SD_182d), " (", num2str(LDD50_SD_182d_lo), ", ", num2str(LDD50_SD_182d_hi), ")")];
        out{4,6} = LDD50_SD_10d/2.7;


        % IT:
        out{6,1} = "IT parameters(kd, mw, Fs, hb)"; out{6,2} = "IT NRMSE"; out{6,3} = "10d LDD50"; out{6,4} = "27d LDD50"; out{6,5} = "182d LDD50"; out{6,6} =  "10d LDD50/2.7"; out{6,7} = "27d LDD10"; out{6,8} = "182d LDD10";
        out{7,1} = par_best_fit_IT([1,2,4,3]); % hb last
        out{7,2} = NRMSE_IT; out{7,3} = strcat(num2str(LDD50_IT_10d), " (", num2str(LDD50_IT_10d_lo), ", ", num2str(LDD50_IT_10d_hi), ")");
        out{7,4} = [strcat(num2str(LDD50_IT_27d), " (", num2str(LDD50_IT_27d_lo), ", ", num2str(LDD50_IT_27d_hi), ")")];
        out{7,5} = [strcat(num2str(LDD50_IT_182d), " (", num2str(LDD50_IT_182d_lo), ", ", num2str(LDD50_IT_182d_hi), ")")];
        out{7,6} = LDD50_IT_10d/2.7;

        %% Print info based on TRT result
        if (SD_TRT_assess || IT_TRT_assess)
            % At least one variant shows TRT (always SD)
            % TRT indicated, use lifespan dose-response. 27d
            % and 182d LDD10s output for sublethal ERA

            out{4,7} = [strcat(num2str(LDD10_SD_27d), " (", num2str(LDD10_SD_27d_lo), ", ", num2str(LDD10_SD_27d_hi), ")")];
            out{4,8} = [strcat(num2str(LDD10_SD_182d), " (", num2str(LDD10_SD_182d_lo), ", ", num2str(LDD10_SD_182d_hi), ")")];

            out{7,7} = [strcat(num2str(LDD10_IT_27d), " (", num2str(LDD10_IT_27d_lo), ", ", num2str(LDD10_IT_27d_hi), ")")];
            out{7,8} = [strcat(num2str(LDD10_IT_182d), " (", num2str(LDD10_IT_182d_lo), ", ", num2str(LDD10_IT_182d_hi), ")")];

            if NRMSE_SD <= NRMSE_IT
                % substance has TRT
                Result = "Step 4: TRT indicated, use GUTS parameters for 27 and 182d effects of all PEQs";
                out{1,1} = "Result: "; out{1,2} = Result;
            else
                % SD shows TRT, but IT fits better
                Result = strcat("Step 3a: Poorer fit variant shows TRT. NRMSE ratio ", num2str(NRMSE_SD/NRMSE_IT));
                out{1,1} = "Result: "; out{1,2} = Result;
            end
        else

            % writecell(out, outTableName);
            Result = "Step 4: No TRT";
            out{1,1} = "Result: "; out{1,2} = Result;
            close all

        end
    end
else
    % Outputs for cases when GUTS cannot be fitted
    if mortFlag
        % Less than 10% mortality. New study required
        Result = "Step 2: TRT assessment cannot be done, insufficient mortality. Specific TRT study required";
        out = cell(1,2);
        out{1} = "Result: "; out{2} = Result;
        close all
    else
        % DRC but for a different reason than earlier
        [LDDs, slope, LD50_emp] = haber_slope_fit(surv, concs, 0.9); % 0.9 for 10% effect

        out = cell(5);
        Result = "Step 2c: Insufficient mortality to fit GUTS, worst case Haber's slope of -2 applied";
        Result2 = "Consider performing a specific TRT study";
        out{1,1} = "Result: "; out{1,2} = Result; out{2,2} = Result2;

        out{4,1} = "DRC slope"; out{4,2} = "DRC 10d LD50 (inflection point)";
        out{4,3} = "27d LDD10"; out{4,4} = "182d LDD10";
        out{5,1} = slope; out{5,2} = LD50_emp;
        out{5,3} = LDDs(1); out{5,4} = LDDs(2);

        close all
    end
end

disp("TRT assessment complete")
% Assessment complete, save results
save_results(path, outFileName.Value, out)


