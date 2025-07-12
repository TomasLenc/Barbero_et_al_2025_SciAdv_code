% Quantify similarity between the RSM obtained from all vs. average
% frontocentral channels.

clear 

%% 

experiment = '2ioi'; 

par = get_par('experiment', experiment); 

subjects = num2cell(par.subjects); 

response = 'eeg'; 

feat_name_resp = 'X'; 
feat_name_stim = 'ratio';

concat_ep = false; 
feat_snr = false; 

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

save_path = fullfile(par.deriv_path, 'summary_comparisons', 'compare_roi_all_front'); 
mkdir(save_path); 


for i_harm_sel=1:length(harm_selections)

    harmonics = harm_selections{i_harm_sel}; 
    harm_str = vec_to_str(harmonics, 'format','%d', 'sep','-'); 

    %% load data for all chan

    roi_name = 'all';         
    avg_chans = false; 


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

    res_all_chan = load(fname); 


    %% load data for front ROI

    roi_name = 'front';         
    avg_chans = true; 


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

    res_front_chan = load(fname); 

    %% save fname 

    fname_save_base = fullfile(save_path, ...
            sprintf('response-%s_feat-%s_featStim-%s_harm-%s_concatEp-%s_snr-%s_roi-allVsFront', ...
                             response, ...
                             feat_name_resp, ...
                             feat_name_stim, ...
                             harm_str, ...
                             jsonencode(concat_ep), ...
                             snr_str ...
                             ));

   %% use individual RDM for tap as a model of RDM for EEG

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
                                          {res_all_chan.res.rdm}, ...
                                          {res_front_chan.res.rdm},  ...
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
                          res_all_chan.res(i_sub).rdm, ...
                          res_front_chan.res(i_sub).rdm, ...
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

    saveas(f, [fname_save_base, '_corrRdmIndRho.svg']); 
    close(f); 




    %% correlation with best model 

    r_all = [res_all_chan.res.r]';
    r_front = [res_front_chan.res.r]'; 

    z_all = r_to_z(r_all); 
    z_front = r_to_z(r_front); 

    c = [135, 175, 201]/255; 

    f = figure('pos', [1517 809 144 193], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select();       

    s = scatter(repmat(0, 1, length(z_front)), z_front, 60, c, 'filled'); 
    hold on ; 
    s = scatter(repmat(1, 1, length(z_all)), z_all, 60, c, 'filled'); 
    plot(repmat([0;1],1,length(z_front)), [z_front, z_all]', '-', 'linew', 2, 'color', c); 
    ax.XLim = [-0.5, 1.5]; 

    ax.XTick = [0, 1]; 
    ax.XTickLabel = {'front', 'all'}; 
    pnl.margin = [10,10,10,10];

    [h, p, ci, stat] = ttest(z_front, z_all);        
    
    bf10 = bf.ttest(z_front, z_all);   % paired samples BF
    
    txt = text(0.5, 1.1*max([z_front; z_all]), ...
               sprintf('r_f_r_o_n_t=%.2f r_a_l_l=%.2f\nBF10=%.2f, \nt(%d)=%.2f, p=%.2g', ...
                    mean(r_front), mean(r_all), bf10, stat.df, stat.tstat, p)); 
    txt.HorizontalAlignment = 'center'; 

    saveas(f, sprintf('%s_ttestRhoFisher.svg', fname_save_base));       

    
    %% distribution of correlations across categorical model boundaries

    r_all_distr = cat(1, res_all_chan.res.r_all_models); 
    r_front_distr = cat(1, res_front_chan.res.r_all_models); 

    z_all_distr = r_to_z(r_all_distr); 
    z_front_distr = r_to_z(r_front_distr); 
    
    r = diag(corr(z_all_distr', z_front_distr', 'type', 'pearson')); 
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

    [h, p, ci, stats] = ttest(z);        

    txt = text(0, max(z)*1.2, sprintf('r = %.2f\np = %.2g', mean(r), p)); 
    txt.HorizontalAlignment = 'center'; 

    ax.XLim = [-0.2, 0.2]; 
    ax.XAxis.Visible = 'off'; 
    ax.YLim = [-0.1, 3.1]; 

    saveas(f, [fname_save_base, '_distrCorr.svg'])

end

% close all



