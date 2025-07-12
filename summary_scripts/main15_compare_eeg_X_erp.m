% Quantify similarity between EEG RSMs obtained from complex and magnitue
% spectra. 

clear 

%% 

experiment = '2ioi'; 

par = get_par('experiment', experiment); 

subjects = num2cell(par.subjects); 

feat_name_stim = 'ratio';

concat_ep = false; 
feat_snr = false; 

roi_name = 'front';         
avg_chans = true; 

harm_selections = {
    [1:6]
}; 

run_stats = true; 

%%

if feat_snr 
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 
else
    snr_str = 'none'; 
end

save_path = fullfile(par.deriv_path, 'summary_comparisons', 'compare_eeg_X_erp'); 
mkdir(save_path); 


for i_harm_sel=1:length(harm_selections)

    harmonics = harm_selections{i_harm_sel}; 
    harm_str = vec_to_str(harmonics, 'format','%d', 'sep','-'); 

    %% load data for X
    
    response = 'eeg'; 
    feat_name_resp = 'X'; 

    subdirs = fullfile(sprintf('response-%s', response), ...
                             sprintf('harm-%s', harm_str), ...
                             sprintf('concatEp-%s_snr-%s', jsonencode(concat_ep), snr_str), ...
                             sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)) ...
                             ); 

    fname = fullfile(par.summary_path, subdirs, ...
            sprintf('response-%s_feat-%s_featStim-%s_harm-%s_concatEp-%s_snr-%s_roi-%s_avgChans-%s_data.mat', ...
                             response, ...
                             feat_name_resp, ...
                             feat_name_stim, ...
                             harm_str, ...
                             jsonencode(concat_ep), ...
                             snr_str, ...
                             roi_name, ...
                             jsonencode(avg_chans) ...
                             )); 
    res_X = load(fname); 

    %% load data for ERP

    response = 'erp'; 
    feat_name_resp = 'X'; 

    subdirs = fullfile(sprintf('response-%s', response), ...
                         sprintf('lfp-%.0fHz', par.filt_erp_cut_low), ...
                         sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)) ...
                         ); 

    fname = fullfile(par.summary_path, subdirs, ...
                     sprintf('response-%s_feat-erp_featStim-%s_lfp-%.0fHz_roi-%s_avgChans-%s_data.mat', ...
                             response, ...
                             feat_name_stim, ...
                             par.filt_erp_cut_low, ...
                             roi_name, ...
                             jsonencode(avg_chans) ...
                             )); 
    res_erp = load(fname); 
    
    %% save fname 

    fname_save_base = fullfile(save_path, ...
            sprintf('response-eeg_featStim-%s_harm-%s_concatEp-%s_snr-%s_roi-%s_avgChans-%s_feat-XvsERP', ...
                             feat_name_stim, ...
                             harm_str, ...
                             jsonencode(concat_ep), ...
                             snr_str, ...
                             roi_name, ...
                             jsonencode(avg_chans) ...
                             ));

    %% correlate RDMs based on X and ERP

    % load stimulus 
    load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 
    rdm_stim = rsa_stim(par, 1/par.cycle_dur * harmonics, ...
                           'feat', feat_name_stim); 

    % group-level permutation test
    if exist([fname_save_base, '_permGroup.mat']) && ~run_stats
        clear res_group
        load([fname_save_base, '_permGroup.mat']); 
    else
        res_group = []; 
        [res_group.r_obs, res_group.r_null_dist, res_group.p_perm] = ...
                                perm_rdm_corr_group_empirical_model(...
                                          {res_X.res.rdm}, ...
                                          {res_erp.res.rdm},  ...
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
                          res_X.res(i_sub).rdm, ...
                          res_erp.res(i_sub).rdm, ...
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

    saveas(f, [fname_save_base, '_corrRdmsIndRho.svg']); 
    close(f); 



    %% correlation with best model 

    r_X = [res_X.res.r]';
    r_erp = [res_erp.res.r]'; 

    z_X = r_to_z(r_X); 
    z_erp = r_to_z(r_erp); 

    c = [135, 175, 201]/255; 

    f = figure('pos', [1517 809 144 193], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select();       

    s = scatter(repmat(0, 1, length(z_X)), z_X, 60, c, 'filled'); 
    hold on ; 
    s = scatter(repmat(1, 1, length(z_erp)), z_erp, 60, c, 'filled'); 
    plot(repmat([0;1],1,length(z_X)), [z_X, z_erp]', '-', 'linew', 2, 'color', c); 
    ax.XLim = [-0.5, 1.5]; 

    ax.XTick = [0, 1]; 
    ax.XTickLabel = {'X', 'ratios'}; 
    pnl.margin = [10,10,10,10];
    
    [h, p, ci, stats] = ttest(z_X, z_erp);        
    
    bf10 = bf.ttest(z_X, z_erp);   % paired samples BF
    
    txt = text(0.5, 1.1*max([z_X; z_erp]), ...
               sprintf('r_X=%.2f r_r_a_t=%.2f\nBF10=%.2f, \np=%.2g', ...
                    mean(r_X), mean(r_erp), bf10, p)); 
    txt.HorizontalAlignment = 'center'; 

    saveas(f, sprintf('%s_ttestRhoFisher.svg', fname_save_base)); 

    %% distribution of correlations across categorical model boundaries

    r_X_distr = cat(1, res_X.res.r_all_models); 
    r_erp_distr = cat(1, res_erp.res.r_all_models); 

    z_X_distr = r_to_z(r_X_distr); 
    z_erp_distr = r_to_z(r_erp_distr); 
    
    r = diag(corr(z_X_distr', z_erp_distr', 'type', 'pearson')); 
    z = r_to_z(r_erp); 

    f = figure('pos', [652 701 85 231], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select(); 

    c = [135, 175, 201] / 255; 
    plot(0 + 0.2*rand(1,length(z))-0.1, z, 'o', 'color', c, 'markerfacecolor', c, 'markersize', 5); 

    z_mu = mean(z); 
    z_sem = std(z) / sqrt(length(z)); 
    z_ci = norminv(1 - 0.025) * z_sem; 

    hold on 
    c = [0 0 0]; 
    plot(0, z_mu, 'o', 'color', c, 'markerfacecolor', c, 'markersize', 8); 
    plot([0,0], [z_mu-z_ci, z_mu+z_ci], 'color', c, 'linew', 3); 

    ax.XLim = [-0.2, 0.2]; 
    ax.XAxis.Visible = 'off'; 

    [h, p, ci, stats] = ttest(z);        

    txt = text(0, max(z)*1.2, sprintf('r = %.2f\np = %.2g', mean(r), p)); 
    txt.HorizontalAlignment = 'center'; 

    ax.YLim = [-0.1, 3.1]; 

    saveas(f, [fname_save_base, '_distrCorr.svg'])

end

close all

