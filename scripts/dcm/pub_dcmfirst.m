%% First level spectral DCM estimation 

%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging
% Pablo Mallaroni, Defined: Katrin Preller, Adeel Razi

%% Initialise
close all
clear all

%% Paths
main_path = '/Users/administrator/Documents/MATLAB/project_cb_tom/dcm';
paths.spm =  '/Users/administrator/Documents/MATLAB/imaging_tools/toolboxes/spm12';
paths.parcs=  '/Users/administrator/Documents/MATLAB/imaging_tools/MNIparcs';
paths.roi = '/Users/administrator/Documents/MATLAB/project_cb_tom/results/tom_rois';

addpath(genpath(main_path),paths.spm, paths.roi,genpath(paths.parcs));

%% Parameters
spm('defaults','fmri'); %128 Hz high pass already set in base defaults
spm_get_defaults('cmdline',true);
spm_jobman('initcfg');

%edit according to your data
subjects = arrayfun(@(x) sprintf('%02d', x), 1:22, 'UniformOutput', false); %sub 1:22
sessions = {'cb','psil','pla'};
%sessions = {'a','b'};
str_func = 'task-rest';

str_roi = 'tom_roi'; 
str_preproc = 'swta';  %string now includes fsltopup after realignment but can be edited
str_dcm = 'dcm_tom_neurosynth_dmn'; 


%% Load coordinate
roi_dat = readtable('/Users/administrator/Documents/MATLAB/project_cb_tom/neurosynth_maps_dmn/neurosynth_coord.csv');
roi_labels = roi_dat.Region;
roi_coords = roi_dat{:,2:end};

%% Begin 
spm
for sub = 1:numel(subjects)
       for ses = 1:numel(sessions)
            % Locate files / directories / aquisition information
            if ~exist(fullfile(main_path,"derivatives",['sub-',subjects{sub}],['ses-',sessions{ses}],'1stlevel_corrected_model_neurosynth',[str_dcm,'.mat']),'file')
            display(['running first-level DCM sub-',subjects{sub},' ses-',sessions{ses}])
            data_path = fullfile(main_path,'nifti_bids',['sub-',subjects{sub}],['ses-',sessions{ses}]); %edit path
            result_path =fullfile(main_path,"derivatives",['sub-',subjects{sub}],['ses-',sessions{ses}]);
            cd(result_path);
            f = spm_select('ExtList',result_path, strcat('^',str_preproc, '.*\.nii$') ,Inf);
            f_info = char(fullfile(data_path, 'func', {dir(fullfile(data_path, 'func', ['*', str_func, '*_bold.json'])).name})); 
            f_dat= fileread(f_info); 
            f_dat= jsondecode(f_dat); %BIDS .JSON aquisition info 

                        
            RT = f_dat.RepetitionTime; %repetition time 
            

            %  Initial GLM for extracting WM/CSF regressors
            %--------------------------------------------------------------------------          
            glmdir = fullfile(result_path,'1stlevel');
            if ~exist(glmdir,'dir')
                mkdir(glmdir)
            end
            
            if ~exist(fullfile(result_path,'1stlevel_corrected_model_neurosynth',[str_dcm,'.mat']),'file')

                clear matlabbatch;
           
                % SPM specification
                matlabbatch{1}.spm.stats.fmri_spec.dir          = cellstr(glmdir);
                matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = RT;
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans   = cellstr(f);
                
                % SPM estimation
                matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glmdir,'SPM.mat'));
                
                % ROI extraction
                matlabbatch{3}.spm.util.voi.spmmat  = cellstr(fullfile(glmdir,'SPM.mat'));
                matlabbatch{3}.spm.util.voi.adjust  = NaN;
                matlabbatch{3}.spm.util.voi.session = 1;
                matlabbatch{3}.spm.util.voi.name    = 'CSF';
                matlabbatch{3}.spm.util.voi.roi{1}.sphere.centre     = [ 0 -40 -5];
                matlabbatch{3}.spm.util.voi.roi{1}.sphere.radius     = 4;
                matlabbatch{3}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
                matlabbatch{3}.spm.util.voi.roi{2}.mask.image        = cellstr(fullfile(glmdir,'mask.nii'));
                matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';
                
                matlabbatch{4} = matlabbatch{3};
                matlabbatch{4}.spm.util.voi.name = 'WM';
                matlabbatch{4}.spm.util.voi.roi{1}.sphere.centre = [0 -30 -25]; 
                
                spm_jobman('run',matlabbatch);
               
                %  Second GLM for extracting ROI timeseries after including WM/CSF regressors
                %--------------------------------------------------------------------------          
                glmdir = fullfile(result_path,'1stlevel_corrected_model_neurosynth');
                if ~exist(glmdir,'dir')
                    mkdir(glmdir); 
                end
                
                clear matlabbatch;
                
                % SPM specification
                matlabbatch{1}.spm.stats.fmri_spec.dir          = cellstr(glmdir);
                matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = RT;
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans   = cellstr(f);
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {
                    spm_select('FPList', result_path,'^rp_av.*\.txt$'),...
                    fullfile(result_path,'1stlevel','VOI_CSF_1.mat'),...
                    fullfile(result_path,'1stlevel','VOI_WM_1.mat'),...
                    }';
                
                % SPM estimation
                matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glmdir,'SPM.mat'));
                
                % ROI extraction
                matlabbatch{3}.spm.util.voi.spmmat  = cellstr(fullfile(glmdir,'SPM.mat'));
                matlabbatch{3}.spm.util.voi.adjust  = NaN;
                matlabbatch{3}.spm.util.voi.session = 1;
                matlabbatch{3}.spm.util.voi.roi{2}.mask.image  = cellstr(fullfile(glmdir,'mask.nii'));
                matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';

                %Loop through parcellation 
                count = 0;
                for roi = 1:numel(roi_labels)
                    matlabbatch{3+count} = matlabbatch{3};
                    matlabbatch{3+count}.spm.util.voi.name = roi_labels{roi};
                    matlabbatch{3+count}.spm.util.voi.roi{1}.sphere.centre = roi_coords(roi,:);
                    matlabbatch{3+count}.spm.util.voi.roi{1}.sphere.radius = 8;
                    matlabbatch{3+count}.spm.util.voi.roi{2}.mask.image    = cellstr(fullfile(glmdir,'mask.nii'));
                    matlabbatch{3+count}.spm.util.voi.expression = 'i1 & i2';
                    count = count+1;
    
                end 

                spm_jobman('run',matlabbatch);
    

                %  SPECIFY & ESTIMATE DCM
                %--------------------------------------------------------------------------         
    
                clear DCM;
                
                % ROIs
                cd(glmdir)
                for roi = 1:numel(roi_labels)
                load(['VOI_',roi_labels{roi},'_1.mat']);
                DCM.xY(roi) = xY; 
                end 
                

                % Metadata
                v = length(DCM.xY(1).u); % number of time points
                n = length(DCM.xY);      % number of regions
                
                DCM.v = v;
                DCM.n = n;
                
                % Timeseries
                DCM.Y.dt  = RT;
                DCM.Y.X0  = DCM.xY(1).X0;
                DCM.Y.Q   = spm_Ce(ones(1,n)*v);
                for i = 1:DCM.n
                    DCM.Y.y(:,i)  = DCM.xY(i).u;
                    DCM.Y.name{i} = DCM.xY(i).name;
                end
    
    
                
                % Task inputs
                DCM.U.u    = zeros(v,1);
                DCM.U.name = {'null'};         
                
                % Connectivity
                DCM.a  = ones(n,n);
                DCM.b  = zeros(n,n,0);
                DCM.c  = zeros(n,0);
                DCM.d  = zeros(n,n,0);
                
                
                % Timing
                DCM.TE     = f_dat.EchoTime;
                DCM.delays = repmat(RT,DCM.n,1);
                
                % Options
                DCM.options.nonlinear  = 0;
                DCM.options.two_state  = 0;
                DCM.options.stochastic = 0;
                DCM.options.analysis   = 'CSD';
                DCM.options.nograph    = 1;
                
                %You can change the name of the DCM here
                str = sprintf(str_dcm);
                DCM.name = str;
                save(fullfile(glmdir,str),'DCM');
                
                DCM = spm_dcm_fmri_csd(fullfile(glmdir,str));
                
            end 
         end
            
       end
end 