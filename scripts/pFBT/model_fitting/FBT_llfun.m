%%%%%%%%%%%%%%%%%
%Giles Story London 2023, adapted from Sam Ereira 2021

% Model function for agent specific learning model

% INPUTS:

% p:    a vector of parameter indices; indexes r.opt_idx to recover
%       full parameter vector
%r:     structure containing input settings, (see FBT_config.m)
%sj:    subject number
%mu:    array of prior means for parameters
%nui:   inverse prior covariance matrix for parameters
%doprior: logical, if set to 1 use prior, if set to 0 do ML

% OUTPUTS:
%
% l:  Posterior probability of data given current parameter setting
% Bs: Model estimated belief about p_Self on each trial
% Bo: Model estimated belief about p_Other on each trial
% RP_hat: Fitted responses, including choice noise

function [l,Bs,Bo,RP_hat]=FBT_llfun(p,r,sj,mu,nui,doprior)

X = r.subjects(sj).data;
Cue = [X().cue];  %trial type - privileged (1), shared (2) or decoy (3)
Obsv = [X().outcome]; %outcome - 0 or 1

%Subject's probability estimates
RP=r.subjects(sj).RP;

%Fitted responses
RP_hat=nan(length(Cue),1);

%Index of probe trials 1=self 2=other
probe=[X().probe];

%Expand parameter vector
P=p(r.p)';

%Transform parameters into bounded space
p_tr =  sigmtr(P,r.LB,r.UB,50);

%Set any params designated as zero to zero
p_tr(r.pfixzero) = 0;

%Name params
alphaS_solo=p_tr(1);
alphaS_shar=p_tr(2);
alphaO_solo=p_tr(3);
alphaO_shar=p_tr(4);
tauS=p_tr(5);
tauO=p_tr(6);
deltaS=p_tr(7);
deltaO=p_tr(8);
lambdaS=p_tr(9);
lambdaO=p_tr(10);

%Calculate prior
if doprior
    x=p;
    lp = -1/2 * (x-mu)'*nui*(x-mu) - 1/2*log(2*pi/det(nui));
    lp=-lp;
else
    lp = 0;
end

%Loop over trials
Scount = 0;
Ocount = 0;
boundary = 0.0001;
Bs(1) = 0.5;
Bo(1) = 0.5;

for n = 1:length(Cue)



    if Cue(n) == 1 %Privileged trial
        PEs(n) = Obsv(n) - Bs(n);
        PEo(n) = 0;
        alphaS = alphaS_solo;
        alphaO = alphaO_solo;
    elseif Cue(n) == 2  %Shared trial
        PEs(n) = Obsv(n) - Bs(n);
        PEo(n) = Obsv(n) -   Bo(n);
        alphaS = alphaS_shar;
        alphaO = alphaO_shar;
    elseif Cue(n) == 3 %Decoy trial
        PEo(n) = Obsv(n) -   Bo(n);
        PEs(n) = 0;
        alphaS = alphaS_solo;
        alphaO = alphaO_solo;
    end

    Bs(n+1) = Bs(n) + alphaS*(PEs(n)+lambdaO*PEo(n)) +  deltaS*(0.5-Bs(n)) ;
    Bo(n+1) = Bo(n) + alphaO*(PEo(n)+ lambdaS*PEs(n)) + deltaO*(0.5-Bo(n)) ;

    %Bound extreme estimates to avoid overflow
    Bs(n+1)=max(Bs(n+1),boundary);
    Bs(n+1)=min(Bs(n+1),1-boundary);

    Bo(n+1)=max(Bo(n+1),boundary);
    Bo(n+1)=min(Bo(n+1),1-boundary);

    %Find likelihood of reported probability from pdf of beta distribution with mode of modelled response
    if probe(n)==1

        Scount = Scount + 1;

        mu_bf = Bs(n+1);

        betaA =  ((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)^3/(27*tauS^3) + (((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)^3/(27*tauS^3) + (- 12*tauS*mu_bf^3 + 16*tauS*mu_bf^2 - 7*tauS*mu_bf + tauS)/(2*tauS) - ((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)*(3*tauS - 14*mu_bf*tauS + 16*mu_bf^2*tauS + mu_bf^2 - 2*mu_bf^3))/(6*tauS^2))^2 - ((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)^2/(9*tauS^2) - (3*tauS - 14*mu_bf*tauS + 16*mu_bf^2*tauS + mu_bf^2 - 2*mu_bf^3)/(3*tauS))^3)^(1/2) + (- 12*tauS*mu_bf^3 + 16*tauS*mu_bf^2 - 7*tauS*mu_bf + tauS)/(2*tauS) - ((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)*(3*tauS - 14*mu_bf*tauS + 16*mu_bf^2*tauS + mu_bf^2 - 2*mu_bf^3))/(6*tauS^2))^(1/3) + ((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)^2/(9*tauS^2) - (3*tauS - 14*mu_bf*tauS + 16*mu_bf^2*tauS + mu_bf^2 - 2*mu_bf^3)/(3*tauS))/((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)^3/(27*tauS^3) + (((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)^3/(27*tauS^3) + (- 12*tauS*mu_bf^3 + 16*tauS*mu_bf^2 - 7*tauS*mu_bf + tauS)/(2*tauS) - ((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)*(3*tauS - 14*mu_bf*tauS + 16*mu_bf^2*tauS + mu_bf^2 - 2*mu_bf^3))/(6*tauS^2))^2 - ((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)^2/(9*tauS^2) - (3*tauS - 14*mu_bf*tauS + 16*mu_bf^2*tauS + mu_bf^2 - 2*mu_bf^3)/(3*tauS))^3)^(1/2) + (- 12*tauS*mu_bf^3 + 16*tauS*mu_bf^2 - 7*tauS*mu_bf + tauS)/(2*tauS) - ((- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)*(3*tauS - 14*mu_bf*tauS + 16*mu_bf^2*tauS + mu_bf^2 - 2*mu_bf^3))/(6*tauS^2))^(1/3) + (- mu_bf^3 + mu_bf^2 - 7*tauS*mu_bf + 3*tauS)/(3*tauS);
        betaB =  (mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)/(3*tauS) + ((mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)^2/(9*tauS^2) - (4*mu_bf + 5*tauS - 18*mu_bf*tauS + 16*mu_bf^2*tauS - 5*mu_bf^2 + 2*mu_bf^3 - 1)/(3*tauS))/((((- 12*tauS*mu_bf^3 + 20*tauS*mu_bf^2 - 11*tauS*mu_bf + 2*tauS)/(2*tauS) - (mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)^3/(27*tauS^3) + ((mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)*(4*mu_bf + 5*tauS - 18*mu_bf*tauS + 16*mu_bf^2*tauS - 5*mu_bf^2 + 2*mu_bf^3 - 1))/(6*tauS^2))^2 - ((mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)^2/(9*tauS^2) - (4*mu_bf + 5*tauS - 18*mu_bf*tauS + 16*mu_bf^2*tauS - 5*mu_bf^2 + 2*mu_bf^3 - 1)/(3*tauS))^3)^(1/2) + (mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)^3/(27*tauS^3) - (- 12*tauS*mu_bf^3 + 20*tauS*mu_bf^2 - 11*tauS*mu_bf + 2*tauS)/(2*tauS) - ((mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)*(4*mu_bf + 5*tauS - 18*mu_bf*tauS + 16*mu_bf^2*tauS - 5*mu_bf^2 + 2*mu_bf^3 - 1))/(6*tauS^2))^(1/3) + ((((- 12*tauS*mu_bf^3 + 20*tauS*mu_bf^2 - 11*tauS*mu_bf + 2*tauS)/(2*tauS) - (mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)^3/(27*tauS^3) + ((mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)*(4*mu_bf + 5*tauS - 18*mu_bf*tauS + 16*mu_bf^2*tauS - 5*mu_bf^2 + 2*mu_bf^3 - 1))/(6*tauS^2))^2 - ((mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)^2/(9*tauS^2) - (4*mu_bf + 5*tauS - 18*mu_bf*tauS + 16*mu_bf^2*tauS - 5*mu_bf^2 + 2*mu_bf^3 - 1)/(3*tauS))^3)^(1/2) + (mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)^3/(27*tauS^3) - (- 12*tauS*mu_bf^3 + 20*tauS*mu_bf^2 - 11*tauS*mu_bf + 2*tauS)/(2*tauS) - ((mu_bf - 4*tauS + 7*mu_bf*tauS - 2*mu_bf^2 + mu_bf^3)*(4*mu_bf + 5*tauS - 18*mu_bf*tauS + 16*mu_bf^2*tauS - 5*mu_bf^2 + 2*mu_bf^3 - 1))/(6*tauS^2))^(1/3);

        SelfPr(Scount) = pdf('beta', RP(n), real(betaA), real(betaB));
        RP_hat(n) = betarnd(real(betaA), real(betaB));

    elseif probe(n)==2

        Ocount = Ocount + 1;

        mu_bf = Bo(n+1);

        betaA =  ((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)^3/(27*tauO^3) + (((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)^3/(27*tauO^3) + (- 12*tauO*mu_bf^3 + 16*tauO*mu_bf^2 - 7*tauO*mu_bf + tauO)/(2*tauO) - ((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)*(3*tauO - 14*mu_bf*tauO + 16*mu_bf^2*tauO + mu_bf^2 - 2*mu_bf^3))/(6*tauO^2))^2 - ((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)^2/(9*tauO^2) - (3*tauO - 14*mu_bf*tauO + 16*mu_bf^2*tauO + mu_bf^2 - 2*mu_bf^3)/(3*tauO))^3)^(1/2) + (- 12*tauO*mu_bf^3 + 16*tauO*mu_bf^2 - 7*tauO*mu_bf + tauO)/(2*tauO) - ((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)*(3*tauO - 14*mu_bf*tauO + 16*mu_bf^2*tauO + mu_bf^2 - 2*mu_bf^3))/(6*tauO^2))^(1/3) + ((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)^2/(9*tauO^2) - (3*tauO - 14*mu_bf*tauO + 16*mu_bf^2*tauO + mu_bf^2 - 2*mu_bf^3)/(3*tauO))/((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)^3/(27*tauO^3) + (((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)^3/(27*tauO^3) + (- 12*tauO*mu_bf^3 + 16*tauO*mu_bf^2 - 7*tauO*mu_bf + tauO)/(2*tauO) - ((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)*(3*tauO - 14*mu_bf*tauO + 16*mu_bf^2*tauO + mu_bf^2 - 2*mu_bf^3))/(6*tauO^2))^2 - ((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)^2/(9*tauO^2) - (3*tauO - 14*mu_bf*tauO + 16*mu_bf^2*tauO + mu_bf^2 - 2*mu_bf^3)/(3*tauO))^3)^(1/2) + (- 12*tauO*mu_bf^3 + 16*tauO*mu_bf^2 - 7*tauO*mu_bf + tauO)/(2*tauO) - ((- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)*(3*tauO - 14*mu_bf*tauO + 16*mu_bf^2*tauO + mu_bf^2 - 2*mu_bf^3))/(6*tauO^2))^(1/3) + (- mu_bf^3 + mu_bf^2 - 7*tauO*mu_bf + 3*tauO)/(3*tauO);
        betaB =  (mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)/(3*tauO) + ((mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)^2/(9*tauO^2) - (4*mu_bf + 5*tauO - 18*mu_bf*tauO + 16*mu_bf^2*tauO - 5*mu_bf^2 + 2*mu_bf^3 - 1)/(3*tauO))/((((- 12*tauO*mu_bf^3 + 20*tauO*mu_bf^2 - 11*tauO*mu_bf + 2*tauO)/(2*tauO) - (mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)^3/(27*tauO^3) + ((mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)*(4*mu_bf + 5*tauO - 18*mu_bf*tauO + 16*mu_bf^2*tauO - 5*mu_bf^2 + 2*mu_bf^3 - 1))/(6*tauO^2))^2 - ((mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)^2/(9*tauO^2) - (4*mu_bf + 5*tauO - 18*mu_bf*tauO + 16*mu_bf^2*tauO - 5*mu_bf^2 + 2*mu_bf^3 - 1)/(3*tauO))^3)^(1/2) + (mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)^3/(27*tauO^3) - (- 12*tauO*mu_bf^3 + 20*tauO*mu_bf^2 - 11*tauO*mu_bf + 2*tauO)/(2*tauO) - ((mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)*(4*mu_bf + 5*tauO - 18*mu_bf*tauO + 16*mu_bf^2*tauO - 5*mu_bf^2 + 2*mu_bf^3 - 1))/(6*tauO^2))^(1/3) + ((((- 12*tauO*mu_bf^3 + 20*tauO*mu_bf^2 - 11*tauO*mu_bf + 2*tauO)/(2*tauO) - (mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)^3/(27*tauO^3) + ((mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)*(4*mu_bf + 5*tauO - 18*mu_bf*tauO + 16*mu_bf^2*tauO - 5*mu_bf^2 + 2*mu_bf^3 - 1))/(6*tauO^2))^2 - ((mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)^2/(9*tauO^2) - (4*mu_bf + 5*tauO - 18*mu_bf*tauO + 16*mu_bf^2*tauO - 5*mu_bf^2 + 2*mu_bf^3 - 1)/(3*tauO))^3)^(1/2) + (mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)^3/(27*tauO^3) - (- 12*tauO*mu_bf^3 + 20*tauO*mu_bf^2 - 11*tauO*mu_bf + 2*tauO)/(2*tauO) - ((mu_bf - 4*tauO + 7*mu_bf*tauO - 2*mu_bf^2 + mu_bf^3)*(4*mu_bf + 5*tauO - 18*mu_bf*tauO + 16*mu_bf^2*tauO - 5*mu_bf^2 + 2*mu_bf^3 - 1))/(6*tauO^2))^(1/3);

        OtherPr(Ocount) = pdf('beta', RP(n), real(betaA), real(betaB)); %find probability of reported probability from pdf of beta with mode of modelled response
        RP_hat(n) = betarnd(real(betaA), real(betaB));
   

    end


end

%Further check to prevent -inf log likelihoods
SelfPr=max(SelfPr,eps);
OtherPr=max(OtherPr,eps);

%Calculate likelihood
ll= min(-sum(log([SelfPr, OtherPr])),10000);  %Upper bound neg log likelihood at 10000

%Posterior
l  = sum(ll) +sum(lp);


