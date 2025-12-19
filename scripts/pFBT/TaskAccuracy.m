%% This script demonstrates how task accuracy was estimated in the pFBT in Mallaroni et al. 2025
%  code written by Sam Ereira (2025)

clear
clc

conditions = {'psilocybin','2cb','placebo'};

data = dir('synthetic_data'); data = {data(3:end).name};

Nsj = 3; %number of synthetic subjects

alpha = 0.15; %arbitrary learning rate for simulating an observer

for icond = 1:3
    for isj = 1:Nsj
     
       load (['synthetic_data' filesep  'synthetic_data_Sj' num2str(isj) '_' conditions{icond} '.mat'])

        %Reports on Self probe trials
        S_reports = [SODdata.trials([SODdata.trials.probe]==1).ReportedProbability];
        
        %Reports on Other probe trials
        O_reports = [SODdata.trials([SODdata.trials.probe]==2).ReportedProbability];
        
        %Ground truth on Self probe trials
        S_truth = [SODdata.trials([SODdata.trials.probe]==1).GroundTruth];
        
        %Ground truth on Other probe trials
        O_truth = [SODdata.trials([SODdata.trials.probe]==2).GroundTruthOther];
        
        %Concatenate response data
        x = [S_reports,O_reports];

        %Loop over sampling trials and implement simplified learning model
        Bs = []; Bs(1) = 0.5;
        Bo = []; Bo(1) = 0.5;
        Scount = 0; Ocount = 0;
        boundary = 0.0001;
        Cue = [SODdata.trials().cue];  %trial type - privileged (1), shared (2) or decoy (3)
        Obsv = [SODdata.trials().outcome]; %outcome - 0 or 1

        for n = 1:length(Cue)
            if Cue(n) == 1 %Privileged trial
                PEs(n) = Obsv(n) - Bs(n); %Self prediction error
                PEo(n) = 0; %Simulated prediction error of the Other
            elseif Cue(n) == 2  %Shared trial
                PEs(n) = Obsv(n) - Bs(n);
                PEo(n) = Obsv(n) -   Bo(n);
            elseif Cue(n) == 3 %Decoy trial
                PEo(n) = Obsv(n) -   Bo(n);
                PEs(n) = 0;
            end

            %Parallel belief updates for Self and Other
            Bs(n+1) = Bs(n) + alpha*(PEs(n));
            Bo(n+1) = Bo(n) + alpha*(PEo(n));

            %Bound extreme estimates to avoid overflow
            Bs(n+1)=max(Bs(n+1),boundary);
            Bs(n+1)=min(Bs(n+1),1-boundary);

            Bo(n+1)=max(Bo(n+1),boundary);
            Bo(n+1)=min(Bo(n+1),1-boundary);
        end

        Bs = Bs(2:end); Bo = Bo(2:end);

        %Only select sampling trials that were immediately followed by
        %probe trials
        Bs = Bs([SODdata.trials.probe]==1);
        Bo = Bo([SODdata.trials.probe]==2);

        %Concatenated simulated choice data for purpose of estimating
        %"optimal" behaviour
        y_sim = [Bs,Bo];

        %Concatenated ground truth data
        y_gt = [S_truth,O_truth];

        %Correlation between real choice data and ground truth
        true_corr_gt = corr(x',y_gt', 'type', 'Spearman');

        %what does chance-level random responding look like
        for inull = 1:1000
            null_corr_gt(inull) = corr(x(randperm(length(x)))', y_gt', 'type', 'Spearman');
        end

        %Accuracy relevative to chance
        balanced_corr_gt(icond,isj) = true_corr_gt - mean(null_corr_gt);

        %"optimal" performance on this task
        optimal_behaviour_corr(icond,isj) = corr(y_sim', y_gt');
    end
end



% Now visualise accuracy
dotsize = 10; dottransparency = 1; error_bar_thick = 2.5; errorcol = [0 0 0];
facecol = [0.6 0.6 0.9]; facecol_psilocybin = [0.6 0.6 0.9];
facecol_2cb =  [0.9 0.6 0.6]; facecol_pla = [0.6 0.8 0.6];
facealpha = 0.3; axisthick =2;

close all
data = balanced_corr_gt;
medians = mean(data)
q1 = std(data)/sqrt(size(data,1));
q3 = std(data)/sqrt(size(data,1));

bar(1,medians(1), 'FaceColor', facecol_psilocybin, 'EdgeColor', 'none', 'FaceAlpha', facealpha);
hold on
bar(2,medians(2), 'FaceColor', facecol_2cb, 'EdgeColor', 'none', 'FaceAlpha', facealpha);
hold on
bar(3,medians(3), 'FaceColor', facecol_pla, 'EdgeColor', 'none', 'FaceAlpha', facealpha);
hold on

% Plot asymmetric error bars (IQR)
hold on;
errorbar(1:size(data,2), medians, q1, q3, 'LineStyle', 'none', 'LineWidth', error_bar_thick, 'CapSize', 0, 'Color', errorcol);
hold off;

for i = 1:size(data,2)
    hold on
    scatter(i+0.25*(rand(size(data,1),1)-0.5), data(:,i), dotsize, 'filled', 'k', 'MarkerFaceAlpha', dottransparency)
end

set(gca, 'FontSize', 14, 'FontName', 'Arial');
ylabel('Performance relative to chance (R)', 'Rotation', 90)

xticks([1:size(data,2)]);                           % Tick positions
xticklabels(conditions);         % Corresponding labels

ax = gca;  
ax.FontWeight = 'bold';

% Set axis line thickness
ax.LineWidth = axisthick;  % Increase the thickness (default is usually 0.5)

box off

low_bound_optimal = prctile(optimal_behaviour_corr(:),75);
upp_bound_optimal = prctile(optimal_behaviour_corr(:),100);


% Define the optimality bounds
y1 = low_bound_optimal; % Lower bound
y2 = upp_bound_optimal; % Upper bound

hold on

% Add dotted horizontal lines without labels
yline(y1, 'k--', 'LineWidth', 0.5);
yline(y2, 'k--', 'LineWidth', 0.5);

% Shade area between the two lines
xLimits = xlim;
fill([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], [y1 y1 y2 y2], ...
     facecol, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Place a left-aligned label inside the shaded area
xText = xLimits(1) + 0.1 * range(xLimits);  % Small offset from left
yText = mean([y1 y2]);  % Center vertically in the shaded region
text(xText, yText, {'Optimal', 'performance'}, ...
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'middle', ...
     'FontWeight', 'normal', ...
     'FontSize', 12)


hold on


ylim([-0.1 0.9])


