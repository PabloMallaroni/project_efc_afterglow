%%%%%%%%%%%%%%%%%
%Giles Story London 2023
%Adapted by Sam Ereira 2024

%Function to configure data and model settings for
% fitting of mixed effects model to probabilistic false belief task data

%Inputs:
%
% model: see EM_fit_SODmod
%options: structure containing options

%Outputs:
%
% r: a structure containing all the necessary settings and input data

function [r] = FBT_config(model, options)

%Get subject ids
directory =  'synthetic_data';
X = dir(directory);
SjList = {X.name}; SjList = SjList(find(contains({X.name}, 'synthetic_data') & contains({X.name}, options.condition )))';

r.nsjs=length(SjList);  %set this as needed
r.SjList = SjList;

%Configure optimiser
try
    r.optr =options.optimiser;
catch
    r.optr ='fminunc';
end
r.options=optimset('display','off','DerivativeCheck','on','MaxFunEvals',10000);

%Set objective function
r.objfunhan = @FBT_llfun;
r.objfun = 'FBT_llfun';

%Configure parameters to optimise and their bounds

%Bounds of all possible parameters

%[alphaS: solo trials,
% alphaS: shared trials
% alphaO: solo trials
% alphaO: shared trials
% tauS
% tauO
% deltaS
% deltaO
% lambdaS
% lambdaO]

%lower and upper bounds of parameters
r.LB = [0 0 0 0 0.0001 0.0001 0 0 -1 -1];
r.UB = [1 1 1 1  0.08  0.08   1 1  1  1];

%Which parameters to optimise
r.model = model;
n_alpha = model(1); %number of learning rate parameters
n_tau = model(2); %number of temperature parmaeters
n_delta = model(3); %number of memory decay paramteres
n_lambda = model(4); %number of leak parameters

%Set the index of params to be optimised - r.opt_idx
pm_ind=1:length(r.UB);
r.opt_idx = [pm_ind(1:n_alpha), pm_ind(5:4+n_tau), ...
    pm_ind(7:6+n_delta) pm_ind(9:8+n_lambda)];

% Reverse parameter index, used to expand out the parameter
% vector again when running the objective function

% r.p lists each param's index in r.opt_idx

if n_alpha == 1 %single learning rate
    p=[1 1 1 1];
elseif n_alpha == 2 %2 agent-specific learning rates
    p=[1 1 2 2];
end
p_tmp = p;
p_tmp(find(p_tmp==999))=[];

if n_tau == 1 %1 temperature parameter
    p = [p max(p_tmp)+1 max(p_tmp)+1];
elseif n_tau == 2 %2 agent-specific temperature parameters
    p= [p max(p_tmp)+1 max(p_tmp)+2];
end
p_tmp = p;
p_tmp(find(p_tmp==999))=[];

if n_delta == 0 %No memory decay
    p= [p 999 999];  %placeholders here, the param will be set to zero in the objective fun
elseif n_delta == 1 %1 shared memory decay parameter
    p = [p max(p_tmp)+1 max(p_tmp)+1];
elseif n_delta == 2 %2 agent-specific memory decay parameters
    p= [p max(p_tmp)+1 max(p_tmp)+2];
end
p_tmp = p;
p_tmp(find(p_tmp==999))=[];


if n_lambda == 0 %No leak parameters
    p=[p 999 999];  %placeholders here, the param will be set to zero in the objective fun
elseif n_lambda == 1  %Symmetric leakage
    p = [p max(p_tmp)+1 max(p_tmp)+1]; %set dlambda to zero
elseif n_lambda == 2 %2 asymmetrical leak parameters
    p= [p max(p_tmp)+1 max(p_tmp)+2];
end

r.p=p;

%Index of which parameters to set to zero
r.pfixzero = find(r.p==999);  %index of any params to fix to zero
r.p(r.p==999)=1; % one is a placeholder, the param will be set to zero in the objective fun

%Specify priors for MAP or first round of EM
init_mu = zeros(1,length(r.LB)); %Prior mean %invsigmtr(r.MB,r.LB,r.UB,50) returns 0
init_sig = 0.25*abs(invsigmtr(0.99*(r.UB-r.LB),r.LB,r.UB,50)); %Prior std

r.init_mu=init_mu(r.opt_idx)';  %note prior means are set to zero in the main code
r.init_sig=init_sig(r.opt_idx)';

init_nu = diag(r.init_sig.^2);  %Prior covariance matrix
init_nui = pinv(diag(r.init_sig.^2));

r.init_nui=init_nui;
r.init_nu=init_nu;

%Set max iterations for EM algorithm
if options.doem
    r.maxit=15;
    if options.doprior_init
        r.fittype='EM_ip';
    else
        r.fittype='EM';
    end
else
    r.maxit=1;
    if options.doprior_init
        r.fittype='MAP';
    else
        r.fittype='ML';
    end
end

%Configure data for each subject
sj=1;
for sji=1:r.nsjs
    ID = SjList{sji};
    load([directory filesep ID]);

    if strcmp(options.fit,'optimum') %Option to fit to ground truth
        sjind=1;
    else
        if strcmp(options.fitsjs,'all')  %Fit to all subjects
            sjind=1;
        elseif strcmp(options.fitsjs,'sel') %Otherwise select sjs based on whether they behave above chance
            sjind=sigsjs(sji);
        end
    end



    X=SODdata.trials;



    if sjind  %If fitting this subject's data
        probe=[X().probe];

        %Extract ground truth on probe trials
        gdtruths=[X.GroundTruth];
        gdtrutho=[X.GroundTruthOther];

        truth=nan(size(gdtruths));
        truth(probe==1)= gdtruths(probe==1);
        truth(probe==2)= gdtrutho(probe==2);

        %RP is 'reported probability' - the response used to fit the model
        if strcmp(options.fit,'optimum') %Option to fit to ground truth - for optimal parameters
            RP=truth';
        elseif strcmp(options.fit,'recovery')  %Option to fit to simulated data for parameter recovery
            RP=options.simsubjects(sj).RP_hat;
        else %Fit to subject's responses
            RP = [X.ReportedProbability]';
        end
        RP=max(RP,eps); %Replace zeros with small number

        %Extra variable for RP and ground truth on probe trials only, used
        %for plotting
        RPin=RP(~isnan(probe));
        truthin=truth(~isnan(probe))';

        %Store in input structure
        r.subjects(sj).data=X;
        r.subjects(sj).probe=probe;
        r.subjects(sj).sjnum=sji;
        r.subjects(sj).RP=RP;
        r.subjects(sj).Nch=sum(~isnan(probe));
        r.subjects(sj).truth=truth;
        r.subjects(sj).gdtruths=gdtruths;
        r.subjects(sj).gdtrutho=gdtrutho;
        r.subjects(sj).truthin=truthin;
        r.subjects(sj).RPin=RPin;

        sj=sj+1;
    end
end

end









