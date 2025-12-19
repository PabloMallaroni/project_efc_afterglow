%%%%%%%%%%%%%%%%%
%Giles Story London 2023

%Function to fit model to probabilistic false belief task data, with
%options for ML, MAP or mixed effects optimisation

% Calls a fitting routine emfit.m

% Inputs
%
% model:    1x4 array specifying which model to fit. The first element
%           specifies the number of learning rates (1 or 2). The second
%           element specifies the number of temperature parameters (1 or 2).
%           The 3rd element specifies the number of memory decay parameters
%           (0, 1 or 2) and the 4th element specifies the number of leak
%           parameters (0, 1 or 2) - if set to 2 there are separate
%           leakage parameters for Self->Other and Other->Self.
%           E.g. EM_fit_FBT([1,1,0,0], options) will fit a model with 1
%           learning rate and 1 temperature parameter.
%
% options:  structure to specify options for fitting
%           options.optimiser - string specifying optimiser to use; if
%           unset defaults to 'fminunc'
%
%           options.fit - string specifying which data to fit to: 'data'
%           fits to subject data. 'optimum' fits to the ground truth to obtain
%           optimal parameters,'recovery' fits to simulated data for parameter recovery
%           (note that 'recovery' requires that options includes a field
%           called simsubjects containing simulated data for each subject)
%
%           options.fitsjs - string specifying which subjects to fit: 'all'
%           fits to all, 'sel' fits to subjects whose performance is above
%           chance.
%
%           options.doem - logical - when set to 1 fits mixed effects (hierarchical)
%           model with empirical priors - set to 0 fits either ML or MAP,
%           depending on the setting of options.doprior_init

%           options.doprior_init - logical - when set to 1 uses a prior on
%           the first round of fitting - if doem set to 0 this selects
%           whether to use MAP (doprior_init=1) or ML (doprior_init=0)

% Outputs

% R:        A structure containing results (see comments below)

function [R] = FBT_fit(model, options)

%Todays date
[d,m,y] = ymd(datetime);
dt = strcat(num2str(y), num2str(m), num2str(d));
ld = [dt '-' num2str(model) ];

%Configure data
[r] = FBT_config(model, options);

%Run models
Np=length(r.opt_idx);
[E,V,~,stats,bf] = FBT_emfit(r,Np,cell(Np,1),2000,0,1,r.maxit,0,'','',options.dostats,options.doprior_init);

%Save in results structure
R.r=r;  %Input data
R.E = E;  %Best fitting parameters for each subject
R.V = V;  %Within subject parameter variance
R.stats = stats;  %Contains further stats for each estimated prior parameter
R.bf = bf;  %Contains estimates of the quantities necessary to compute Bayes  factors for model comparison; R.bf.iL is the negative log model evidence

switch options.fit
    case 'optimum'
        save([ld '_' r.fittype '_' options.condition '-opt'],'R');
    case 'recovery'
        save([ld '_' r.fittype  '_' options.condition '-rec'],'R');
    case 'data'
        save([ld '_' r.fittype  '_' options.condition],'R');
end

return






