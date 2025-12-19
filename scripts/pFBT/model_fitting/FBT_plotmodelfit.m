%Giles Story London 2023

%Script to plot simple summary of model fits -
%plots results from the first subject

Figure1=figure;

%%% Data from one subject
sj=12;
%Load model parameters
r=R.r;
L_self =  sigmtr(R.E(4,:),-1,1,50)';
L_other =  sigmtr(R.E(5,:),-1,1,50)';

%Extract reported probabilities on probe trials
X=r.subjects(sj);
probe=[X().probe];
RPs=X.RP(probe==1);
RPo=X.RP(probe==2);

%Extract fitted beliefs
eval(['[~,Bs,Bo]=' r.objfun '(r.p,r,sj,0,0,0);']);
Bs=Bs(probe==1);
Bo=Bo(probe==2);

%Index of self and other probe trials
probetrs=probe(probe==1|probe==2);
s_tr=find(probetrs==1);
o_tr=find(probetrs==2);

%Plot
h=subplot(1,1,1,'Parent',Figure1);
cla(h)
plot(h,s_tr,Bs,'--','Color',[0.909 0.364709 0.5450],'LineWidth',1); hold on
plot(h,s_tr,RPs,'Color',[0.909 0.364709 0.5450],'LineWidth',2); hold on
plot(h,o_tr,Bo,'--','Color',[0.30196 0.745 0.933],'LineWidth',1); hold on
plot(h,o_tr,RPo,'Color',[0.30196 0.745 0.933],'LineWidth',2); hold on
xlabel(h,'\fontsize{12} Probe Trial')
ylabel(h,'\fontsize{12} Estimated P')
title(h,['\fontsize{18} \fontsize{14} \rm Subject ' num2str(sj) ' \lambda_S_e_l_f=' sprintf('%.2g',L_self(sj)) ' \lambda_O_t_h_e_r=' sprintf('%.2g',L_other(sj))])
legend('Model Self','Self','Model Other','Other')









