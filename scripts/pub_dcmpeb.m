%% Second level within-condition spectral DCM estimation 

%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging
% Pablo Mallaroni, Samuel Eriera,  Defined: Katrin Preller, Adeel Razi
%% Initialise
close all; clear all;

main_path     = '/Users/administrator/Documents/MATLAB/project_cb_tom';
paths.spm     = '/Users/administrator/Documents/MATLAB/imaging_tools/toolboxes/spm-main/spm-main';
paths.functions = '/Users/administrator/Documents/MATLAB/imaging_tools/functions';

addpath(genpath(main_path));
addpath(paths.spm);
addpath(paths.functions);

spm('defaults','fmri');
spm_get_defaults('cmdline',true);
spm_jobman('initcfg');

sessions    = {'cb','psil','pla'};
contrasts   = [1;1;0];  % session contrast (drug effect)
str_dcm     = 'dcm_tom_neurosynth_dmn';
regionNames = {'vmPFC','dmPFC','PCu','left TPJ','right TPJ'};
numRegions  = numel(regionNames);
factorNames = {'condition effect','lambda effect'};  % intercept, Y

%% Load behavioural  data (They match - cross check bad but on purpose)
clear Y subjects
dat_table = readtable(fullfile(main_path,'data','behaviour','lambda.xlsx'));
Y = dat_table.leak_diff;
bad_sub = find(isnan(Y));
Y(bad_sub,:) = [];

% Build subjects from the filesystem
d = dir(fullfile(main_path, 'dcm', 'derivatives', 'sub-*'));
d = d([d.isdir]);
subs = {d.name}';

% Sort numerically
subnums = str2double(erase(subs, 'sub-'));
[~, order] = sort(subnums);
subs    = subs(order);
subnums = subnums(order);

% Exclude only 15 and 17
exclude = [15, 17];
subjects = subs(~ismember(subnums, exclude));


%% Within subject second level peb
% Within-subject second-level PEB (spectral DCM)

spm
PEBs = cell(numel(subjects),1);

for s = 1:numel(subjects)
    GCM = cell(numel(sessions),1);
    for ses = 1:numel(sessions)
        result_path = fullfile(main_path,'dcm','derivatives', subjects{s}, ...
                               ['ses-',sessions{ses}], '1stlevel_corrected_model_neurosynth');
        dcm_file = fullfile(result_path, [str_dcm '.mat']);
        if ~exist(dcm_file,'file'); error('Missing DCM file: %s', dcm_file); end
        S = load(dcm_file, 'DCM');
        GCM{ses,1} = S.DCM;
    end

    X = [ones(numel(contrasts),1), contrasts(:)];
    M = struct('Q','single','X',X);
    PEBs{s,1} = spm_dcm_peb(GCM, M, {'A'});
end


%% % Now run "third" level

%create a 3rd level design matrix for 3rd level, where each row is a different subject
X3 = [];
X3(:,1) = ones(length(PEBs),1);
X3(:,2) = Y; %This is where the behavioural parameter of interest goes as a column vector (i.e. self/other leak drug-pla effect)

% %We can add in further columns to X3 as well if we want to control for
% %other covariates, e.g. age, sex, session order
X3(:,2:end) = (X3(:,2:end)- mean(X3(:,2:end))); %demean


M = struct();
M.Q      = 'single';
%M.maxit  = 512; 
M.X      = X3;

%fit PEB (2nd level DCM model with all connections)
%The 3rd level PEB has (N^2 x 2) x P effective connectivity parameters (where N
%is the number of ROIs, and P is the number of columns in the design matrix
%X3). E.g. If we have 4 ROIs and 3 columns in design matrix X3 (a column of
%1s, self-other leakage drug-pla effect, and framewise displacement)
%then we end up with 32 x 3 effective connectivity parameters
%Within each column of 32 EC parameters, the first 16 are probably less interesting, whilst the
%2nd set of 16 pertain to the drug-induced effect and these are of interest. 
%
% Parameters 1-16 in column 1 will show the average effective connectivity
% across all sessions, across subjects
%
% Parameters 17-32 in column 1 will show the average drug-induced change in effective
% connectivity, across subjects 
% 
% Parameters 1-16 in column 2 will show how the average effective
% connectivity across all sessions, varies as a function of the drug-induced behavioural effect. 
% If the BMA retains any of these parameters one way of
% interpreting this is that there are some baseline trait effective
% connectivity phenotypes that predispose people to a stronger behavioural
% drug effect 
%
% Parameters 17-32 in column 2 will show how the drug-induced change in
% effective connectivity varies as a function of the drug-induced
% behavioural effect. 
%
% Parameters in column 3 and beyond will not be of interest. 

PEB3 = spm_dcm_peb(PEBs,M,{'A'}); %Here is where we could also select specific connections (A1,3) etc) 'A(1,1)'
BMA3 = spm_dcm_peb_bmc(PEB3);

%Use this to review BMA results
spm_dcm_peb_review(BMA3); 

out_model_path = fullfile(main_path,'dcm');
%save(fullfile(out_model_path,[str_dcm,'_within_results.mat']),"PEBs","PEB3","BMA3")


%% Plot results.
numFactors = size(X3,2);
eFC  = BMA3.Ep;
pProb =  BMA3.Pp;


clear eFC_matrices pProb_matrices
[eFC_matrices,pProb_matrices] = deal({})
for i = 1:numFactors
    baseIdx = (i-1) * numRegions^2 * 2;  % Start index for this factor
    
    % Reshape eFC and pProb for "all" and "drug" conditions
    eFC_matrices{i, 1} = reshape(full(eFC(baseIdx + 1 : baseIdx + numRegions^2)), numRegions, numRegions);
    eFC_matrices{i, 2} = reshape(full(eFC(baseIdx + numRegions^2 + 1 : baseIdx + 2 * numRegions^2)), numRegions, numRegions);
    
    pProb_matrices{i, 1} = reshape(full(pProb(baseIdx + 1 : baseIdx + numRegions^2)), numRegions, numRegions);
    pProb_matrices{i, 2} = reshape(full(pProb(baseIdx + numRegions^2 + 1 : baseIdx + 2 * numRegions^2)), numRegions, numRegions);
end

% --- make the figure (before the loop)
fig = figure('Color','w','Position',[100 100 1800 1000]);

climsPerRow = [0.3, 4,4];
for i = 1:numFactors
    L = climsPerRow(i);

    subplot(numFactors,2,2*i-1);
    plot_eFC_mat_supplement(eFC_matrices{i,1}, pProb_matrices{i,1}, regionNames, coolwarm(256), L);
    title(['average eFC factor: ' factorNames{i}]);

    subplot(numFactors,2,2*i);
    plot_eFC_mat_supplement(eFC_matrices{i,2}, pProb_matrices{i,2}, regionNames, coolwarm(256), L);
    title(['drug-modulated eFC ' factorNames{i}]);
end

% --- save as SVG (vector)
out_dir = '/Users/administrator/Documents/MATLAB/project_cb_tom/data/behaviour/pub_outputs';
if ~exist(out_dir,'dir'), mkdir(out_dir); end
out_svg = fullfile(out_dir, 'eFC_factor_panels.svg');
%print(fig, out_svg, '-dsvg', '-painters');      % vector SVG

% ----Plot a change matrix we can affix with the same dimensions anyways
% and then threshold retrospectvely
fig = figure('Color','w','Position',[100 100 1800 1000]);
subplot(2,2,1)
plot_eFC_delta(eFC_matrices{1,1}, pProb_matrices{1,1}, eFC_matrices{1,2}, pProb_matrices{1,2}, regionNames, coolwarm(256), 0.3)
% --- save as SVG (vector)
out_dir = '/Users/administrator/Documents/MATLAB/project_cb_tom/data/behaviour/pub_outputs';
out_svg = fullfile(out_dir, 'eFC_factor_panels_diff.svg');
%print(fig, out_svg, '-dsvg', '-painters');      % vector SVG


