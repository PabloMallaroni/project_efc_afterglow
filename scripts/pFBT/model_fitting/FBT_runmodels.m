
%Giles Story London 2023, adapted by Sam Ereira 2024

%Script to fit a set of models to probabilistic false belief task data, with
%options for ML, MAP or mixed effects optimisation

%Set up options as described in FBT_fit.m
clear
clc

conditions = {'placebo', 'psilocybin', '2cb'}
cond = 2; %1 for placebo, 2 for psilocybin, 3 for 2CB - %choose a condition

options.fit='data';
options.doem=0; %EM (expectation-maximisation algorithm is NOT used)
options.doprior_init=1; %MAP (maximum a-posteriori) IS used
options.fitsjs='all';
options.session=1;
options.dostats=1;
options.recovery = 0;
options.condition = conditions{cond};

%Run models with various combinations of parameters

for alpha= [1,2] %1 shared or 2 independent learning rate parameters
    for tau = [2] %1 shared or 2 independent temperature parameters
        for delta= [2] %0, 1 shared or 2 independent memory decay parameters
            for lambda=[2]% 0, 1 shared or 2 independent leak parameters
                model=[alpha tau,delta lambda];
                [R] = FBT_fit(model, options);
            end
        end
    end
end


%FBT_plotmodelfit

