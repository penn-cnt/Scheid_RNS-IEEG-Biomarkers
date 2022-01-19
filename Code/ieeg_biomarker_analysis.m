%% IEEG Biomarker Analysis Pipeline

% Scheid BH, Bernabei JM, Khambhati AN, Mouchtaris S, Jeschke J,
% Bassett DS, et al. Intracranial electroencephalographic biomarker 
% predicts effective responsive neurostimulation for epilepsy prior to 
% treatment. Epilepsia. 2022 Jan 7; DOI: 10.1111/epi.17163
%
% Author: Brittay Scheid (bscheid@seas.upenn.edu)
% December 2021

% ** Instructions ** %
%
% Gather Data:
%  * Create an account on ieeg.org
%  * Download the ieeg cli (https://bitbucket.org/ieeg/ieeg/wiki/Downloads#IEEG_CLI.md)
%  * Add 'ieeg' to your path as a link to ieeg.bin (usually added as /usr/bin/ieeg.bin) 
%  * Use the following commands to access the data
%       * ieeg download-rec-objects --dataset UCSF_Biomarkers_Data UCSF_shared.mat
%       * ieeg download-rec-objects --dataset NYU_Biomarkers_Data NYU_shared.mat
%       * ieeg download-rec-objects --dataset Penn_Biomarkers_Data Penn_shared.mat
%
% Set parameters: 

p.fig_pth = 'Figures';
p.data_pth = '/path/to/shared_data/';

load(fullfile(p.data_pth,'UCSF_shared.mat'))
load(fullfile(p.data_pth,'Penn_shared.mat'))
load(fullfile(p.data_pth,'NYU_shared.mat'))

addpath('Functions')

allSync = [sync_UCSF_shared, sync_Penn_shared, sync_NYU_shared];
allPts = [patients_UCSF_shared, patients_Penn_shared, patients_NYU_shared];

cols = [[0, 118, 192]; [103, 119, 26]; [235, 103, 33]; [163, 2, 52]; [255, 223, 118]]./255; 


%% ANALYSIS SETTINGS 

out_yrs = [1, 2, 3];                    % years since stimOndate to get outcomes
targlabel = sprintfc('%.02f', out_yrs); % label of outcome dates
outTypes =  {'closest'};                %{'closest', 'locf', 'interpolate'};
respThresh = 50;                        % Responders defined as sz percent reduction >= respThresh
dozscore =  false;                      % set to true if pt curves should be zscored
lim= [15, 83];                          % limits for outcome bins

note = '_raw_post_onset';

metricStruct = allSync;             % global network metric structure 
metricName = 'synchronizability';   % metric name   

dispOn = true;    
saveOn = true;     % save figures (needs display to be true)
warning('off', 'stats:perfcurve:SubSampleWithMissingClasses')

% Calculation of score: mean of synchyronizability values after seizure onset 
scoreCalc = (@(x) mean(x(:,round(end/2*1.1):end),2, 'omitNan')); % Take mean of onset

bandlabels=[ {'beta'}, {'highgamma'},  {'broadband_CC'}];
bandlabelnames = [{'\beta band'}, {'high-\gamma band'}, {'broad band'}];

[~,~, all_outcomes] = getOutcomeFromTimepoint([allPts.stimOnDate], {allPts.outcome}, 10);
static_outcomes = cellfun(@(x) all(x>=respThresh)|all(x<respThresh), all_outcomes);

% list of patient populations of interest
patient_sets = {
    {allPts([1:length(allPts)]), 'all patients'}; ...                       % All patients 
    {allPts(static_outcomes), 'static outcomes'};...                        %static outcomes
    {allPts(~static_outcomes), 'dynamic outcomes'};...                      %dynamic outcomes 
    {allPts(contains({allPts.ID}, 'HUP')), 'HUP'}; ...                      % HUP
    {allPts(contains({allPts.ID}, 'NY')), 'NYU'};...                        % NYU only
    {allPts(contains({allPts.ID}, 'NP')), 'UCSF'}; ...                      % UCSF only
    {allPts(contains([allPts.leadLocations], 'N')), 'Neocortical'};...      % Neocortical only
    {allPts(contains([allPts.leadLocations], 'M')),'mesial temporal'};...   % Mesial Temporal only
    {allPts(~contains([allPts.leadLocations], 'N')), '1+ mesial leads'};... % Neocortical only
    {allPts(~contains([allPts.leadLocations], 'M')),'1+ neocortical leads'};... % Neocortical only
    {allPts(contains([allPts.laterality], 'B')), 'bilateral'};              % Bilateral
    {allPts(~contains([allPts.laterality], 'B')), 'lateral'};               % lateral
    {allPts(arrayfun(@(y)sum(strcmp(string(cellfun(@(x)sort(upper(x(1))),...
        allPts(y).IEEG_elecType)), 'D'))./length(allPts(y).IEEG_elecType),...
        (1:30))>=.5), 'majority depth'};...                                  % Majority Depth
    {allPts(arrayfun(@(y)sum(strcmp(string(cellfun(@(x)sort(upper(x(1))),...
        allPts(y).IEEG_elecType)), 'D'))./length(allPts(y).IEEG_elecType),...
        [1:30])<.5), 'majority grid and strip'};...                          $ Majority Grid/Strip
    {allPts(arrayfun(@(y)sum(strcmp(string(cellfun(@(x)sort(upper(x(1))),...
        allPts(y).IEEG_elecType)), 'G'))./length(allPts(y).IEEG_elecType),...
        [1:30])>.5), 'grid only'}};                                           % Grid Only

%% Run all Patient Sets 
rng(30);
mxy = 0; miny = Inf;
for ptset = 1:length(patient_sets)

    select_pts = patient_sets{ptset}{1}; 
    dset = patient_sets{ptset}{2};
    
    figsuf1 = sprintf('%s_%s_zscore-%s%s', dset, metricName, string(dozscore), note);
      
    for outt = outTypes(1)                % Run all types of outTypes
        outType = outt{1};
    
    % Get Outcome labels for Select Pts
    [outcomes, dates, all_outcomes, all_dates] = getOutcomeFromTimepoint([select_pts.stimOnDate],...
    {select_pts.outcome}, out_yrs, 'type', outType); 
    
    AUCs = nan(length(out_yrs), length(bandlabels), 3);
    ptable = nan(length(out_yrs), length(bandlabels),2); 
    mn_diffs = nan(length(out_yrs), length(bandlabels),2); % mean differences and std
    
    %mean(x(round(end/2*1.5):end))
    all_pt_mets = cell2mat(cellfun(@(x) scoreCalc(x), {metricStruct.(metricName)}, 'Uni' , 0))'; 
    [metricStruct.outcome] = deal(nan(size(outcomes,2),1));

if dispOn
    figure(1); clf; tiledlayout('flow'); set(figure(1), 'Position', [65   142   413   739]);  
    figure(2); clf; tiledlayout('flow'); set(figure(2), 'Position', [479   142   464   739])
    figure(3); clf; tiledlayout('flow'); set(figure(3), 'Position', [944   142   390   739]); 
    figure(22); clf; tiledlayout('flow'); set(figure(3), 'Position', [10  142   390   739]); 
end

for i_targ = 1:length(out_yrs)      
    
    if dispOn
        figure(1); axx(i_targ) = nexttile; hold on;
        figure(2); nexttile(); hold on;
        figure(3); ax(i_targ) = nexttile; hold on;
        figure(22); ax(i_targ) = nexttile; hold on;
    end
    
    pt_mets = struct(); ctr = 1;
    for band= 1:length(bandlabels) 
        band_inds_metstruct = strcmp([metricStruct.band], bandlabels{band});
        
        for i_pt = 1:length(select_pts)
            
            pt_inds = strcmp(select_pts(i_pt).ID, [metricStruct.ID]); 
            inds = find(pt_inds.*band_inds_metstruct);
    
            pt_mets(ctr).ID = select_pts(i_pt).ID;
            pt_mets(ctr).band = bandlabels(band);
            pt_mets(ctr).score = nanmean(all_pt_mets(inds));
            pt_mets(ctr).outcome = outcomes(i_pt,:)';
            
            [metricStruct(inds).outcome] = deal(outcomes(i_pt,:)');
                        
            ctr = ctr + 1; 
        end

        i_resp = [pt_mets.outcome] >= 50;
        band_inds =  find(strcmp([pt_mets.band], bandlabels{band}));
        i_resp = i_resp(i_targ, band_inds);

        resp = [pt_mets(band_inds(i_resp)).score];
        nresp = [pt_mets(band_inds(~i_resp)).score];
        labels = [true(1,length(resp)),false(1,length(nresp))];

        [pval, v, stats] = ranksum(resp, nresp); scr = stats.ranksum;
        ptable(i_targ, band,:) = [pval, scr]; 
        
         % Flip AUC curve for case where responders have higher sync value than non-responders
        [X,Y,~,AUC] = perfcurve(labels, [resp, nresp], 0, 'NBoot',500, 'XVals','All');
        if AUC(1) < 0.5
            [X,Y,~,AUC] = perfcurve(labels, [resp, nresp], 1, 'NBoot',500, 'XVals','All');
        end
        AUCs(i_targ, band,:) = AUC; 
               
        mn_diffs(i_targ, band,1) = mean(nresp) - mean(resp);
        mn_diffs(i_targ, band, 2) = norm(std(nresp)/sqrt(length(nresp)),std(resp)/sqrt(length(resp))); % diff b/w sample means

        if dispOn
            figure(1); 
            g1 = band - .2;
            scatter((g1) * ones(sum(i_resp),1), resp, 30, 'markerfacecolor', cols(1,:), 'markeredgecolor', 'black', 'jitter', 'On', 'jitterAmount', .08)
            plot([-.1,.1]+ g1, ones(1,2)*nanmean(resp), 'color', 'black', 'linewidth', 1)
            plot([g1, g1], nanmean(resp)+[-nanstd(resp),nanstd(resp)] , 'color', 'black', 'linewidth', 1)
            g2 = band + .2;
            scatter((g2) * ones(sum(~i_resp),1), nresp,'markerfacecolor', cols(4,:), 'markeredgecolor', cols(4,:), 'jitter', 'On', 'jitterAmount', .08)
            plot([-.1,.1]+ g2, ones(1,2)*nanmean(nresp), 'color', 'black', 'linewidth', 1)
            plot([g2, g2], nanmean(nresp)+[-nanstd(nresp),nanstd(nresp)] , 'color', 'black', 'linewidth', 1)
           
            if pval < 0.05
                plot(band, max([resp,nresp])+.1, '*', 'color', 'black')
                text(band+.05, max([resp,nresp])+.1, sprintf('p = %0.3f', pval)) 
                plot([band-.3, band+.3], ones(1,2)*max([resp,nresp])+.05, 'color', 'black', 'linewidth', 1)
            end

            figure(2); hold on; % Regression Image
            outs = [pt_mets.outcome]; 
            yy = outs(i_targ, band_inds); xx = [pt_mets(band_inds).score];
            [rcorr,pcorr]=corr(xx', yy');
             S= polyfit(xx, yy, 1);
             scatter(xx,yy, 50, cols(band,:), 'filled')
             xs=get(gca, 'xlim'); 
             plot(xs, S(2)+xs*S(1), 'Color', cols(band,:))
             num = randi(length(yy)-band);
             text(xx(band+num)*1.1, yy(band+num)*.9,sprintf('r=%0.2f, p = %0.2f',...
                 rcorr, pcorr), 'fontSize', 13)
             
             figure(22); hold on; % Regression Image box plots... 
             grp = [double(1*outs(i_targ, band_inds)>=lim(2)) + ...
                 2*(outs(i_targ, band_inds)<lim(2) & outs(i_targ, band_inds)> lim(1)) + ...
                 3 * (outs(i_targ, band_inds)<=lim(1))] + (band - 1) * (length(bandlabels) + 1);
             grpMn = groupsummary(xx', grp','mean');
             boxchart(grp, xx, 'BoxFaceColor', cols(band,:), 'WhiskerLineColor',...
                 cols(band,:), 'MarkerColor', cols(band,:))
             plot(unique(grp), grpMn, '-o', 'color', 'black', 'linewidth', 1, 'HandleVisibility', 'off')
             tabulate(grp)
             
             figure(3); hold on; %AUCs
             plot(X,Y(:,1), 'color', cols(band,:));
             auc_resp_plot = [Y(:,2); flip(Y(:,3))];
             %fill([X; flip(X)], auc_resp_plot , 1,'facecolor',cols(band, :),'edgecolor','none', 'facealpha', 0.3, 'HandleVisibility','off');
        end
    end
    
    if dispOn
        figure(1); axis tight
        title(sprintf('%0.2f years', out_yrs(i_targ)))
        ylabel(metricName)
        xticks([1:length(bandlabelnames)]); xticklabels(bandlabelnames)
        xlim([.5, length(bandlabels)+.5])
        mxy = max([mxy, get(gca, 'YLim')]);
        miny = min([miny, get(gca, 'YLim')]);
        
       figure(2); title(sprintf('%0.2f years', out_yrs(i_targ)))
       
       figure(3);
       plot([0,1], [0,1], 'black--','HandleVisibility','off');
       title(sprintf('%s ROC, %s, %s', metricName, dset, targlabel{i_targ}))
       legend(strcat(bandlabelnames, sprintfc('- %.02f', AUCs(i_targ,:,1))), 'Location', 'Best')    
       
       figure(22); 
       title(sprintf('Outcome Regression, %0.2f', out_yrs(i_targ)))
       xticks([1:4*length(bandlabels)-1])
       xticklabels(repmat({sprintf('≥%d%%',lim(2)), sprintf('%d-%d%%', [lim(2), lim(1)]),...
           sprintf('≤%d%%',lim(1)), ''}, 1, length(bandlabels)));
       xtickangle(45)
       xlabel('percent seizure reduction'); ylabel('synchronizability')
       legend(bandlabelnames,'Location', 'Best')
       
    end  
end
    end
    

figure(4); clf; subplot(1,2,1); hold on;
imagesc(squeeze(ptable(:,:,1))); set(gca,'YDir','normal'); colorbar
xticks([1:4]); xticklabels(bandlabels); xlabel('band')
yticks([1:length(out_yrs)]); yticklabels(out_yrs); ylabel('out yrs');
[xx,yy] = find(squeeze(ptable(:,:,1))<=0.05); if xx, plot(yy,xx, '*'), end; axis tight;
title('pvalues')

subplot(1,2,2); 
heatmap(squeeze(AUCs(:,:,1)));
title('AUC')

figure(5); clf; hold on % Change over time
errorbar(repmat(out_yrs',1,3), squeeze(mn_diffs(:,:,1)), squeeze(mn_diffs(:,:,2)),...
    'HandleVisibility', 'Off')
set(gca, 'ColorOrder', cols(1:3,:))
plot(out_yrs, squeeze(mn_diffs(:,:,1)), '-o');
title(sprintf('Mean difference in %s over time', metricName))
ylabel(sprintf('\\Delta %s', metricName)); 
%xticks([1:length(out_yrs)]);
xlabel('years since stimulation on')
legend(bandlabelnames, 'location', 'best')
xlim([min(out_yrs)-.2, max(out_yrs)+.2])   

if saveOn
    
    if ~exist(fullfile(p.fig_pth, metricName), 'dir')
        mkdir(fullfile(p.fig_pth, metricName))
    end
    
    figure(1); 
    for i_targ = 1:length(out_yrs)
        set(axx(i_targ),'Ylim', [miny-0.01, mxy+0.01])
    end   
    set(figure(1), 'Position', [65   142   413   739])
    saveas(figure(1), fullfile(p.fig_pth, metricName, sprintf('scatterPlot_%s.svg', figsuf1)), 'svg')
    saveas(figure(1), fullfile(p.fig_pth, metricName, sprintf('scatterPlot_%s.fig', figsuf1)))
    
    set(figure(2), 'Position', [479   142   464   739])
    saveas(figure(1), fullfile(p.fig_pth, metricName, sprintf('outcomeRegression_%s.svg', figsuf1)), 'svg')
    saveas(figure(1), fullfile(p.fig_pth, metricName, sprintf('outcomeRegression_%s.fig', figsuf1)))
    
    set(figure(22), 'Position', [479   123   347   758])
    saveas(figure(22), fullfile(p.fig_pth, metricName, sprintf('boxplot_outcome_%s.svg', figsuf1)), 'svg')
    saveas(figure(22), fullfile(p.fig_pth, metricName, sprintf('boxplot_outcome_%s.fig', figsuf1)))
    
    figure(3); 
    for i_targ = 1:length(out_yrs)
        legend(ax(i_targ),'Location', 'Best')
    end
    saveas(gcf, fullfile(p.fig_pth, metricName, sprintf('ROC_%s.svg', figsuf1)), 'svg')
    saveas(gcf, fullfile(p.fig_pth, metricName, sprintf('ROC_%s.fig', figsuf1)))

    set(figure(5), 'Position', [1201, 66, 411 , 277])
    saveas(figure(5), fullfile(p.fig_pth, metricName, sprintf('diffOverTime_%s.svg', figsuf1)), 'svg')
    saveas(figure(5), fullfile(p.fig_pth, metricName, sprintf('diffOverTime_%s.fig', figsuf1)))
    
end


%if saveOn
    ps_table = array2table(reshape(ptable, [], length(bandlabels)*2),...
        'variablenames', [strcat(bandlabels, '-p'), strcat(bandlabels, '-U')],...
        'rownames', sprintfc('%0.2f yrs', out_yrs));
    auc_table = array2table(reshape(AUCs, [], length(bandlabels)*3), ...
        'variablenames', [strcat('AUC-', bandlabels), strcat('AUC-low', bandlabels), strcat('AUC-high', bandlabels)], ...
        'rownames', sprintfc('%0.2f yrs', out_yrs));
    
    writetable(ps_table,fullfile(p.fig_pth,'raw_met_results.xls'),'Sheet',dset,'Range','A1', 'WriteRowNames',true)
    writetable(auc_table,fullfile(p.fig_pth,'raw_met_results.xls'),'Sheet',dset,'Range', ...
        sprintf('A%d', length(out_yrs)+3), 'WriteRowNames',true) 
%end


end

%%
figure(3); clf; hold on;
band = 1;
band_inds_metstruct = strcmp([metricStruct.band], bandlabels{band});

ctr = 1;
for i_sz = find(band_inds_metstruct)
    plot(metricStruct(i_sz).synchronizability(1:round(end/2))+ctr)
    ctr = ctr + 1;
end

    