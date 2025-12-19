function [E,V,alpha,stats,bf,fitparams] = FBT_emfit(r,Np,varargin)
%
% %Giles Story Max Planck UCL Centre London 2020 - adapted code from Q Huys
% and D Schad- adds option to use priors for the first round of EM or fit as
% ML or MAP - edits are marked GWS
%
% Note that in Mallaroni et al. 2025, MAP is used (NOT EM)
%
% [E,V,alpha,stats,bf,fitparams] = EMFIT(llfunc,D,Np,[reg],[Nsample],[docheckgrad],[nograd],[maxit],[dofull],[savestr],[loadstr]);
%
% Perform a random-effects fit using expectation-maximimization.
%
% NOTE: This is in development. NO GUARANTEES for correctness...
%
% LLFUNC	- this is a string that points to a likelihood function of choices (a
% model). The function must have the form:
%
%    [l,dl] = llfunc(x,D,musj,nui,doprior);
%
% where x are the parameters, D(sj) is the data for subject sj. MUSJ(:,sj) is
% the prior mean parameter for subject sj. If no regressors (see below) are
% included, then this should be the same for all subjects. If a regressor is
% included for a GLM analysis, then each subject's value will depend on this.
% DOPRIOR determines whether a Gaussian prior with mean musj and diagonal
% covariance matrix NU is added to the log likelihood of the choices.  The
% function must output the sum of the total log likelihood of the choice (l),
% and the gradient wrt all the paramters (dl). For example llfunc see llrw.m
% included at the end of this function
%
% D contains all the data, with D(sj) containing all the data for subject
% sj=1...Nsj. If D(sj).Nch are the number of observations for subjects sj, then
% bf can return integrated BIC values and integrated Laplacian values (see
% below).
%
% NP is th enumber of parameters in the function llfunc.
%
% REG is a cell structure that must be the length NP.  For each parameter, REG
% can be defined. For instance, for a model with two parameters, one might want
% to asj whether the second parameter is related to some psychometric variable
% psi at the group level, taking into account how noisy each subject is. To do
% this, define
%
%     reg{1} = [];
%     reg{2} = psi;
%
% with psi(i) being the psychometric score for subject i.
%
% PARALLEL: If a matlabpool has been opened prior to calling the function it
% will use it (using PARFOR) and hence run faster. Ideally, you should use the
% same number of processors as you have separate subjects in your dataset, or at
% least around 80%.
%
% The OUTPUT is as follows:
%     E        is a matrix of size NpxNsj containing the MAP-EM parameters
%     V        is a matrix of size NpxNsj containing the inverse Hessians around
%              individual subject parameter estimates
%     alpha    contains the estimated coefficients (both means and regression
%              coefficients if REG has been included)
%     stats    contains further stats for each estimated prior parameter, in particular
%              standard error estimates (from finite differences), t and p
%              values, and ML estimates as stats.EML. stats.EX is 1 if it
%              converged; 0 if it reached the MAXIT limit, -2 if a saddle point
%              (with some positive group Hessian diagonals) rather than a maximum was
%              reached, and -1 otherwise.
%     bf       contains estimates of the quantities necessary to compute Bayes
%    			   factors for model comparison. The integrated, group-level BIC
%    			   bf.ibic simply counts the parameters at the group level, while
%    			   bf.ilap uses the finite difference Hessian around the prior mean
%    			   for a Laplacian estimate, which is typically less conservative.
%     fitparam If asjed for, this outputs various parameters used by emfit for
%              fitting - including the data (this might take up a lot of space!)
%
% Additional variables can be defined:
%
%     NSAMPLES is the number of samples used for integration (default: 2000).
%     DOCHECKGRAD=1 sets a flag to check the gradients dl provided by llfunc
%     (default 0). NOGRAD has to be set to 1 if no gradients are provided by
%     llfunc. MAXIT is the maximum number of EM iterations (default 500). If this
%     is set to 0, then ML and only a single EM iteration will be computed. The
%     ML parameters are always in stats.EML. If DOFULL is set to zero, only the
%     diagonal of the Hessians are used for EM. This effectively imposes a prior
%     that sets off-diagonal elements to zero, but is no longer recommended.
%     If SAVESTR is provided, intermediate results are saved in this file.
%
% Copyright Quentin Huys and Daniel Schad, 2015
% www.quentinhuys.com/code.html
% www.quentinhuys.com/pub.html
% qhuys@cantab.net


%=====================================================================================
% setting up
rng('default'); %GWS
dx= 0.001; 													% step for finite differences
fitparams.version='0.151110';							% version of this script

nargin = length(varargin);
t=varargin;
if nargin>0 & ~isempty(t{1}); reg         = t{1}; else reg=cell(Np,1);   end;
if nargin>1 & ~isempty(t{2}); Nsample     = t{2}; else Nsample = 2000;   end;
if nargin>2 & ~isempty(t{3}); docheckgrad = t{3}; else docheckgrad = 0 ; end;
if nargin>3 & ~isempty(t{4}); nograd      = t{4}; else nograd = 0 ;      end;
if nargin>4 & ~isempty(t{5}); maxit       = t{5}; else maxit = 500;      end;
if nargin>5 & ~isempty(t{6}); dofull      = t{6}; else dofull = 1;       end;
if nargin>6 & ~isempty(t{7}); savestr     = t{7}; else savestr = '';     end;
if nargin>7 & ~isempty(t{8}); loadstr     = t{8}; else loadstr = '';     end;
if nargin>8 & ~isempty(t{9}); dostats     = t{9}; else dostats=1;        end; %GWS
if nargin>9 & ~isempty(t{10}); doprior_init= t{10}; else doprior_init=1;        end; %GWS
% Deal with gradients being provided or not
if nograd													% assume gradients are supplied
    fminopt=optimset('display','off','Algorithm','active-set');
else
    fminopt=optimset('display','off','GradObj','on');
    if docheckgrad; 											% if they are, then can check them.
        fminopt=optimset('display','off','GradObj','on','DerivativeCheck','on');
    end
end
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

% check if regressors have been provided correctly
if Np~=length(reg);
    error('You must provide a regressor cell entry for each parameter.');
elseif ~iscell(reg);
    error('REG must be a cell strucure of length Np.');
end

fstr=str2func(r.objfun);									% GWS prepare function string
Nsj = length(r.subjects); sjind=1:Nsj;							% number of subjects
Xreg = repmat(eye(Np),[1 1 Nsj]);
Nreg=0;
for j=1:Np
    if size(reg{j},2)==Nsj; reg{j} = reg{j}';end
    for k=1:size(reg{j},2)
        Nreg = Nreg+1;
        Xreg(j,Np+Nreg,:) = reg{j}(:,k);
    end
end

Npall= Np+Nreg;
coeff_vec = Inf*ones(Npall,1);


%=====================================================================================
fprintf('\nStarting EM estimation');

alpha = zeros(Npall,1); 								% individual subject means
for sj=1:Nsj
    musj(:,sj) = Xreg(:,:,sj)*alpha;
end
%nui = 0.01*eye(Np); nu = inv(nui);					% prior variance over all params
nui = r.init_nui; nu = r.init_nu;					% GWS prior variance over all params
E = zeros(Np,Nsj); % individual subject parameter estimates
emit=0;nextbreak=0;stats.ex=-1;PLold= -Inf;

% continue previous fit
if ~isempty(loadstr);
    eval(['load ' loadstr ' E V alpha stats emit musj nui']);
end
lastPL = 1e9;

if maxit>1
    niters=1;  %if doing EM
else
    niters=10;  %if doing ML or MAP
end

while 1;emit=emit+1;
    % E step...........................................................................
    t0=tic;
    % if docheckgrad; checkgrad(func2str(fstr),randn(Np,1),.001,r(1),musj(:,1),nui,doprior), end
    parfor_progress(Nsj);
    parfor sj=sjind %parfor
        %disp(sj);
        tt0=tic; ex=-1;tmp=0;
        est=[]; hess=[]; fval=[];
        while ex<1           
            %GWS - if doing only one iteration, i.e. MAP or ML fitting
            % if dopriorinit=1 this does MAP;if dopriorinit=0 does ML

            if emit==1  %On the first round use a prespecified prior to prevent wild ML estimates
                doprior=doprior_init;                
            else
                doprior=1;           
            end
        
            for iter=1:niters
                init=r.init_mu +  randn(Np,1).*r.init_sig; 
                [est(iter,:),fval(iter,:),ex(iter,:),~,~, hess(:,:,iter)] = fminunc(@(x)fstr(x,r,sj,musj(:,sj),nui,doprior),init,r.options);
            end

            [fval, min_ind] = min(fval);   %find the minimum over all runs
            est = est(min_ind,:);
            ex = ex(min_ind,:);
            hess = hess(:,:,min_ind);            
        
            if ex<0 ; tmp=tmp+1; fprintf('didn''t converge %i times exit status %i\r',tmp,ex); end

            try
                pinv(full(hess));
            catch
                ex = -1;
            end
        end
        E(:,sj) 		= est;								% Subjets' parameter estimate
        W(:,:,sj) 	= pinv(full(hess));				% covariance matrix around parameter estimate
        V(:,sj) 		= diag(W(:,:,sj));				% diagonal undertainty around parameter estimate
        PL(sj)  		= fval;								% posterior likelihood
        tt(sj) 		= toc(tt0);
        parfor_progress;
        %    fprintf('Emit=%i subject %i exit status=%i\r',emit,sj,ex)
    end
    parfor_progress(0);
    % wait(1000);
    Estore(emit,:,:) = E;
    
    if sum(PL) > lastPL
        %keyboard
    end
    lastPL = sum(PL);
    fprintf('\nEmit=%i PL=%.2f full loop=%.2gs one subj=%.2gs parfor speedup=%.2g',emit,sum(PL),toc(t0),mean(tt),mean(tt)*Nsj/toc(t0))
    if emit==1; stats.EML = E; end
    if emit==2; stats.EMAP0 = E; end
    if nextbreak==1; break;end
    
    % M step...........................................................................
    if emit> 1; 							% only update priors after first MAP iteration
        while 1								% iterate over prior mean and covariance updates until convergence
            % prior mean update - note prior mean is different for each subject if
            % GLM is included
            alpha_old = alpha;
            ah = zeros(Npall,1);
            vh = zeros(Npall);
            for  sj=1:Nsj
                ah = ah + Xreg(:,:,sj)'*nui*E(:,sj);
                vh = vh + Xreg(:,:,sj)'*nui*Xreg(:,:,sj);
            end
            alpha = pinv(vh)*ah; %identical to mean(E,2) if GLM is excluded
            for sj=1:Nsj
                musj(:,sj) = Xreg(:,:,sj)*alpha;
            end
            alpha_store(emit,:) = alpha;
            figure(1); plot(alpha_store); legend;
            figure(2); for pm=1:size(E,1); subplot(2,6,pm);hist(E(pm,:)); end
            % prior covariance update
            if ~dofull 											% use diagonal prior variance
                nunew = diag(sum(E.^2 + V - 2*E.*musj + musj.^2,2)/(Nsj-1));
            else													% use full prior covariance matrix
                nunew = zeros(Np);
                for sj=sjind;
                    nunew = nunew + E(:,sj)*E(:,sj)' - E(:,sj)*musj(:,sj)' ...
                        -musj(:,sj)*E(:,sj)' + musj(:,sj)*musj(:,sj)' + W(:,:,sj);
                end
                nunew = nunew/(Nsj-1);
            end
            if det(nunew)<0; fprintf('negative determinant - nufull not updated');
            else           ; nu = nunew;
            end
            nui = pinv(nu);
            
            if norm(alpha-alpha_old)<1e6; break;end
        end
    end
    
    % check for convergence of EM procedure or stop if only want ML / MAP0..............
    if maxit==1 | (maxit==2 & emit==2); break; end		% stop if only want ML or ML& MAP0
    par(emit,:) = mean(musj,2);		
   
    if emit>1
        critpar=sum(abs(par(emit,:)-par(emit-1,:)))
        if abs(sum(PL)-PLold)<(0.0025*Nsj)||critpar<5
            nextbreak=1;
            stats.ex=1;
            fprintf('...converged');
        end
    end  %GWS - changed to relax convergence criterion
  
    if emit==maxit; nextbreak=1;stats.ex=0;fprintf('...maximum number of EM iterations reached');end
    PLold=sum(PL);
    if length(savestr)>0; eval(['save ' savestr ' E V alpha stats emit musj nui']);end
end
stats.PL = PL;
stats.subjectmeans= musj;
stats.groupvar= nu;

if dostats
    %=====================================================================================
    fprintf('\nComputing individual subject BIC values');
    for sj=sjind
        stats.LL(sj)  = fstr(E(:,sj),r,sj,musj(:,sj),nui,0);
        try
            bf.bic(sj) =  -2*stats.LL(sj)   + Np*log(r.subjects(sj).Nch); 
            bf.aic(sj) =  -2*stats.LL(sj)   + Np*2;
        catch
            bf.bic(sj) =  NaN;
        end
    end
    
   % if maxit<=2; return; end								% end here if only want ML or ML & MAP0
    
    %=====================================================================================
    fprintf('\n');
    
    oo = ones(1,Nsample);
    LLi=zeros(Nsample,1);
    Hess = zeros(Npall,Npall,Nsj);
    for sj=sjind;
        fprintf('\rSampling subject %i ',sj)
        
        % estimate p(choices | prior parameters) by integrating over individual parameters
        muo = musj(:,sj)*oo;
        es = sqrtm(nu)*randn(Np,Nsample)+muo;
        parfor k=1:Nsample; %parfor
            LLi(k) = fstr(es(:,k),r,sj,musj(:,sj),nui,0);
        end
        lpk0 = max(-LLi); 
        pk0 = exp(-LLi-lpk0);
        bf.iL(sj) = log(mean(exp(-LLi-lpk0)))+lpk0;	% integrated likelihood
        
        bf.SampleProbRatio(sj)=(sum(-LLi)/-stats.PL(sj))/Nsample;
        bf.EffNsample(sj)=sum(exp(-LLi/r.subjects(sj).Nch));	% eff # samples
        
        % shift samples to get gradients and (full) Hessian around prior parameters
        des = es-muo;
        err = sum(des.*(nui*des),1);
        for l=1:Npall
            % shift *samples* by +2*dx
            foo = alpha; foo(l) = foo(l)+2*dx; mud = Xreg(:,:,sj)*foo;
            desshift = es - mud*oo;
            lw = -1/2*sum(desshift.*(nui*desshift),1) - (-1/2*err);
            w = exp(lw'); w = w/sum(w);
            lldm = log(pk0'*w)+lpk0;
            
            % shift *samples* by +dx
            foo = alpha; foo(l) = foo(l)+dx; mud = Xreg(:,:,sj)*foo;
            desshift = es - mud*oo;
            lw = -1/2*sum(desshift.*(nui*desshift),1) - (-1/2*err);
            w = exp(lw'); w = w/sum(w);
            lld(l,1) = log(pk0'*w)+lpk0;
            bf.EffNsampleHess(l,l,sj) = 1/sum(w.^2);
            
            % finite difference diagonal Hessian elements
            Hess(l,l,sj) = (lldm+bf.iL(sj)-2*lld(l,1))/dx^2;
            
            % now compute off-diagonal terms
            for ll=1:l-1
                % again shift samples, but along two dimensions
                foo  = alpha;  foo(l) = foo(l)+dx; foo(ll) = foo(ll)+dx; mud = Xreg(:,:,sj)*foo;
                desshift = es - mud*oo;
                lw = -1/2*sum(desshift.*(nui*desshift),1) - (-1/2*err);
                w = exp(lw'); w = w/sum(w);
                lldd(l,ll) = log(pk0'*w)+lpk0;  % off-diagonal second differential
                bf.EffNsampleHess(l,ll,sj) = 1/sum(w.^2);	% eff # samples
                bf.EffNsampleHess(ll,l,sj) = 1/sum(w.^2);
                
                % off-diagonal Hessian terms
                Hess(l,ll,sj) = (lldd(l,ll) - lld(l,1) - lld(ll,1) + bf.iL(sj))/dx^2;
                Hess(ll,l,sj) = Hess(l,ll,sj);
            end
        end
        if any(any(bf.EffNsampleHess(:,:,sj)<50))
            warning('Warning: Less than 50 effective samples - dimensionality prob.  too high!');
        end
    end
    fprintf('...done ')
    
    stats.individualhessians = Hess;
    stats.groupmeancovariance = - pinv(sum(Hess,3))*Nsj/(Nsj-1);
    saddlepoint=0;
    if any(diag(stats.groupmeancovariance)<0);
        warning('Negative Hessian, i.e. not at maximum - try running again, increase MAXIT if limit reached')
        stats.ex=-2;
        saddlepoint=1;
    end
    
    % compute t and p values comparing each parameter to zero
    if ~saddlepoint												% can estimate p values
        stats.groupmeanerr	= sqrt(diag(stats.groupmeancovariance));
        stats.tval				= alpha./stats.groupmeanerr;
        stats.p					= 2*tcdf(-abs(stats.tval), Nsj-Npall);
    else																% can't estimate p values
        stats.groupmeanerr	= NaN*alpha;
        stats.tval				= NaN*alpha;
        stats.p					= NaN*alpha;
    end
    
    %=====================================================================================
    fprintf('\nComputing iBIC and iLAP')
    
    Nch=0; for sj=sjind; Nch = Nch + r.subjects(sj).Nch;end
    bf.ibic =  -2*(sum(bf.iL) - 1/2*(2*Np+Nreg)*log(Nch));
    bf.ilap =  -2*(sum(bf.iL) - 1/2*   Np      *log(Nch)  + 1/2*log(det(stats.groupmeancovariance)));
    
    if nargout>=6
        %=====================================================================================
        fprintf('\nSaving fit parameters ')
        fitparams.likelihoodfunction=llfunc;
        fitparams.reg=reg;
        fitparams.Nsample=Nsample;
        fitparams.docheckgrad=docheckgrad;
        fitparams.nograd=nograd;
        fitparams.maxit=maxit;
        fitparams.dofull=dofull;
        fitparams.r=r;
        fitparams.Np=Np;
    end
    fprintf('\nDone\n ')
    
    return
else
    bf=[];
    stats=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


