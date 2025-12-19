%% Extract parameters from winning model

clear
clc
close all

load 'Parameter Estimates.mat'
%Data struct R contains raw parameter estimates and model fit statistics
%for the Mallaroni 2025 dataset
%vectors alpha, lambda and tau contain the raw native-space parameter estimates for the N =
%20 subjects included in the final anlysis

%transform alpha from native space into model space (bounded by 0 and 1)
alpha = sigmtr(alpha,0,1,50);

%transform tau from native space into model space (bounded by 0.001 and 0.08)
tau = sigmtr(tau,0.001,0.08,50);

%transform lambda from native space into model space (bounded by -1 and 1)
lambda = sigmtr(lambda,-1,1,50);




%% Linear mixed effect model for statistical inference

IDs = 1:length(alpha) %subject ID

leak = abs(lambda).*alpha; %self-other mergence factor
leak_placebo = leak(:,3);
leak_2cb = leak(:,2);
leak_psi = leak(:,1);

n = size(leak,1); % number of participants

% Vectors of all data concatenated
leak_all = [leak_placebo; leak_2cb; leak_psi];

% Subject IDs repeated for each condition
subject = repmat((1:n)', 3, 1);

% Condition labels
condition = [repmat("placebo", n, 1);
             repmat("2cb",   n, 1);
             repmat("psi",   n, 1)];

% Create table
T = table(subject, condition, leak_all, ...
          'VariableNames', {'Subject', 'Condition', 'Leak'});

% Create a new logical variable: 1 = drug, 0 = placebo
T.IsDrug = ismember(T.Condition, ["2cb", "psi"]);
T.Condition = categorical(T.Condition, {'placebo', '2cb', 'psi'});
T.IsDrug = categorical(T.IsDrug); % so fitlme can handle it properly

lme = fitlme(T, 'Leak ~ IsDrug + (Condition|Subject)')

stats = anova(lme)


%% Visualise Self-Other mergence (leak) factor

dotsize = 10; dottransparency = 1; error_bar_thick = 2.5; errorcol = [0 0 0];
facecol = [0.6 0.6 0.9]; facecol_psilocybin = [0.6 0.6 0.9];
facecol_2cb =  [0.9 0.6 0.6]; facecol_pla = [0.6 0.8 0.6];
facealpha = 0.3; axisthick =2;

tmp = leak; 
medians = nanmedian(tmp)
q1 = nanstd(tmp)/sqrt(size(tmp,1));
q3 = nanstd(tmp)/sqrt(size(tmp,1));

bar(1,medians(1), 'FaceColor', facecol_psilocybin, 'EdgeColor', 'none', 'FaceAlpha', facealpha);
hold on
bar(2,medians(2), 'FaceColor', facecol_2cb, 'EdgeColor', 'none', 'FaceAlpha', facealpha);
hold on
bar(3,medians(3), 'FaceColor', facecol_pla, 'EdgeColor', 'none', 'FaceAlpha', facealpha);
hold on

% Plot asymmetric error bars (IQR)
hold on;
errorbar(1:size(medians,2), medians, q1, q3, 'LineStyle', 'none', 'LineWidth', error_bar_thick, 'CapSize', 0, 'Color', errorcol);
hold off;

for i = 1:size(tmp,2)
    hold on
    scatter(i+0.25*(rand(size(tmp,1),1)-0.5), tmp(:,i), dotsize, 'filled', 'k', 'MarkerFaceAlpha', dottransparency)
end

set(gca, 'FontSize', 14, 'FontName', 'Arial');
%yticks(0:0.02:0.08); 
ylabel('Self-Other Mergence', 'Rotation', 90)

xticks([1 2 3]);                           % Tick positions
xticklabels({'Psilocybin', '2-CB', 'Placebo'});         % Corresponding labels

ax = gca;  
ax.FontWeight = 'bold';

% Set axis line thickness
ax.LineWidth = axisthick;  % Increase the thickness (default is usually 0.5)

box off
