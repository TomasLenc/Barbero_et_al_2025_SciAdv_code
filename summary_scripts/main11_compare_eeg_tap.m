% Quantify similarity between the EEG and tapping RSM of each subject. 

clear 

response_tap = 'tap-force'; 

feat_name_eeg = 'X'; 
feat_name_tap = 'X'; 
feat_name_stim = 'ratio';

roi_name = 'all';         
avg_chans = false; 

% roi_name = 'front';         
% avg_chans = true; 

concat_ep = false; 
feat_snr = false; 

harm_selections = {
    {[1:6], [1:12]}
}; 

run_stats = true; 


%%

experiment = '2ioi'; 

par = get_par('experiment', experiment); 

subjects = num2cell(par.subjects); 

%%

if feat_snr 
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 
else
    snr_str = 'none'; 
end


for i_harm_sel=1:length(harm_selections)

    harmonics_eeg = harm_selections{i_harm_sel}{1}; 
    harmonics_tap = harm_selections{i_harm_sel}{2}; 

    harm_str_eeg = vec_to_str(harmonics_eeg, 'format','%d', 'sep','-'); 
    harm_str_tap = vec_to_str(harmonics_tap, 'format','%d', 'sep','-'); 

    %% load eeg

    subdirs = fullfile(sprintf('response-%s', 'eeg'), ...
                             sprintf('harm-%s', harm_str_eeg), ...
                             sprintf('concatEp-%s_snr-%s', jsonencode(concat_ep), snr_str), ...
                             sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)) ...
                             ); 

    fname_load_base_eeg = fullfile(par.summary_path, subdirs, ...
            sprintf('response-%s_feat-%s_featStim-%s_harm-%s_concatEp-%s_snr-%s_roi-%s_avgChans-%s', ...
                             'eeg', ...
                             feat_name_eeg, ...
                             feat_name_stim, ...
                             harm_str_eeg, ...
                             jsonencode(concat_ep), ...
                             snr_str, ...
                             roi_name, ...
                             jsonencode(avg_chans) ...
                             )); 

    res_eeg = load([fname_load_base_eeg, '_data.mat']); % res

    %% load tap
    
    save_path = fullfile(par.deriv_path, 'summary_comparisons', 'compare_eeg_tap'); 

    if strcmp(feat_name_tap, 'X')
        subdirs = fullfile(sprintf('response-%s', response_tap), ...
                                 sprintf('harm-%s', harm_str_tap), ...
                                 sprintf('concatEp-%s_snr-%s', jsonencode(concat_ep), snr_str) ...
                                 ); 

        fname_load_base_tap = fullfile(par.summary_path, subdirs, ...
                sprintf('response-%s_feat-%s_featStim-%s_harm-%s_concatEp-%s_snr-%s', ...
                                 response_tap, ...
                                 feat_name_tap, ...
                                 feat_name_stim, ...
                                 harm_str_tap, ...
                                 jsonencode(concat_ep), ...
                                 snr_str)); 

        save_path = fullfile(save_path, sprintf('response-eeg-%s', response_tap), ...
                                 sprintf('concatEp-%s_snr-%s', jsonencode(concat_ep), snr_str), ...
                                 sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)) ...
                                 ); 

    elseif strcmp(feat_name_tap, 'ratio')

        subdirs = fullfile(sprintf('response-%s', 'tap-onsets')); 

        fname_load_base_tap = fullfile(par.summary_path, subdirs, ...
                                       sprintf('response-tap-onsets'));  

        save_path = fullfile(save_path, sprintf('response-eeg-tap-onsets'), ...
                             sprintf('concatEp-%s_snr-%s', jsonencode(concat_ep), snr_str), ...
                             sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)) ...
                             ); 
    else
        error('bad tap feature name')
    end

    res_tap = load([fname_load_base_tap, '_data.mat']); % res

    %% prepare save-path 

    mkdir(save_path); 

    fname_save_base = fullfile(save_path, ...
            sprintf('featEEG-%s_featTAP-%s_featStim-%s', ...
                             feat_name_eeg, ...
                             feat_name_tap, ...
                             feat_name_stim)); 
    if strcmp(feat_name_eeg, 'X')
        fname_save_base = sprintf('%s_harmEEG-%s', fname_save_base, harm_str_eeg); 
    end
    if strcmp(feat_name_tap, 'X')
        fname_save_base = sprintf('%s_harmTAP-%s', fname_save_base, harm_str_tap); 
    end
    if strcmp(feat_name_stim, 'X')
        fname_save_base = sprintf('%s_harmSTIM-%s', fname_save_base, harm_str_eeg); 
    end

    %% use individual RDM for tap as a model of RDM for EEG

    % load stimulus 
    load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 
    rdm_stim = rsa_stim(par, 1/par.cycle_dur * harmonics_eeg, ...
                           'feat', feat_name_stim); 

    % group-level permutation test
    if exist([fname_save_base, '_permGroup.mat']) && ~run_stats
        clear res_group
        load([fname_save_base, '_permGroup.mat']); 
    else
        res_group = []; 
        [res_group.r_obs, res_group.r_null_dist, res_group.p_perm] = ...
                                perm_rdm_corr_group_empirical_model(...
                                          {res_eeg.res.rdm}, ...
                                          {res_tap.res.rdm},  ...
                                          par.n_perm_rsa_group, ...
                                          'parfor', true, ...
                                          'covariate', rdm_stim); 
        save([fname_save_base, '_permGroup.mat'], 'res_group'); 
    end        

    % individual level permutation test with fixed model
    if exist([fname_save_base, '_permInd.mat']) && ~run_stats
        clear res_ind
        load([fname_save_base, '_permInd.mat']); 
    else
        res_ind = []; 
        for i_sub=1:length(subjects)
            i_sub
            [res_ind(i_sub).r_obs, res_ind(i_sub).p_perm, res_ind(i_sub).r_null_dist] = ...
                perm_rdm_corr(...
                          res_eeg.res(i_sub).rdm, ...
                          res_tap.res(i_sub).rdm, ...
                          par.n_perm_rsa_group, ...
                          'parfor', true, ...
                          'covariate', rdm_stim); 
        end
        save([fname_save_base, '_permInd.mat'], 'res_ind'); 
    end

    % RHO of individual subjects 
    f = figure('pos', [809 696 80 173], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select(); 

    plot_r_bar([res_ind.r_obs], [res_ind.p_perm], 'y_maxlim', 1, 'ax', ax); 

    ax.YLim = [-1, 1]; 
    ax.YTick = [-1, 1]; 
    ax.YAxis.LineWidth = 2; 
    ax.XAxis.LineWidth = 2; 
    txt = text(0, 1, sprintf('p=%.1g', res_group.p_perm)); 
    txt.HorizontalAlignment = 'center'; 
    pnl.margin = [8,4,0,1]; 

    saveas(f, [fname_save_base, '_corrRdmEegTapIndRho.svg']); 
    close(f); 


    %% distribution of correlations across categorical model boundaries
    
    r_eeg = cat(1, res_eeg.res.r_all_models); 
    r_tap = cat(1, res_tap.res.r_all_models); 

    z_eeg = r_to_z(r_eeg); 
    z_tap = r_to_z(r_tap); 

    r = diag(corr(z_eeg', z_tap', 'type', 'pearson')); 
    z = r_to_z(r); 

    z_mu = mean(z); 
    z_sem = std(z) / sqrt(length(z)); 
    z_ci = norminv(1 - 0.025) * z_sem; 

    f = figure('pos', [652 701 85 231], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select(); 
    c = [135, 175, 201] / 255; 
    plot(0 + 0.2*rand(1,length(z))-0.1, z, 'o', 'color', c, 'markerfacecolor', c, 'markersize', 5); 
    hold on 
    c = [0 0 0]; 
    plot(0, z_mu, 'o', 'color', c, 'markerfacecolor', c, 'markersize', 8); 
    plot([0,0], [z_mu-z_ci, z_mu+z_ci], 'color', c, 'linew', 3); 
    ax.XLim = [-0.2, 0.2]; 
    ax.XTick = []; 
    ax.XAxisLocation = 'origin'; 
    ax.YAxis.LineWidth = 2; 
    ax.XAxis.LineWidth = 2; 
    ax.XAxis.TickLength = [0,0]; 
    ax.YAxis.TickLength = [0,0]; 
    ax.YTick = sort([min(ax.YTick), 0, max(ax.YTick)]); 

    [h, p, ci, stats] = ttest(z);        

    txt = text(0, max(z)*1.2, sprintf('r = %.2f\np = %.2g', mean(r), p)); 
    txt.HorizontalAlignment = 'center'; 

    ylim([-inf, max(z)*1.2])
    ax.YLim = [ceil(min(min(z), 0)*100)/100, floor(max(z)*100)/100];         
    ax.YTick = [ceil(min(min(z), 0)*100)/100, floor(max(z)*100)/100];         

    saveas(f, [fname_save_base, '_eegTapDistrCorr.svg'])
    close(f); 

end



% close all


