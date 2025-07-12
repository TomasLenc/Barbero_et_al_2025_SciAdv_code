% Summarize RSA results for the auditory nerve model. 

clear 

%% 

experiment = '2ioi'; 

par = get_par('experiment', experiment); 

feat_name_resp = 'X'; 
feat_name_stim = 'ratio';

response = 'urear-an'; 

concat_ep = false; 
feat_snr = false; 

harm_selections = {
    [1:6]
    [1:12]
}; 

run_stats = false; 

plot_big_figure = false; 

%%

if feat_snr 
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 
else
    snr_str = 'none'; 
end

for i_harm_sel=1:length(harm_selections)

    harmonics = harm_selections{i_harm_sel}; 

    %% prepare output paths

    harm_str = vec_to_str(harmonics, 'format','%d', 'sep','-'); 

    subdirs = fullfile(sprintf('response-%s', response), ...
                             sprintf('harm-%s', harm_str), ...
                             sprintf('concatEp-%s_snr-%s', jsonencode(concat_ep), snr_str) ...
                             ); 

    save_path = fullfile(par.summary_path, subdirs); 
    mkdir(save_path); 

    fname_save_base = fullfile(save_path, ...
            sprintf('response-%s_feat-%s_featStim-%s_harm-%s_concatEp-%s_snr-%s', ...
                             response, ...
                             feat_name_resp, ...
                             feat_name_stim, ...
                             harm_str, ...
                             jsonencode(concat_ep), ...
                             snr_str ...
                             )); 

    %% load data 

    fname_tbl = sprintf('response-%s_feat-%s_featStim-%s_harm-%s_concatEp-%s_snr-%s.csv', ...
                             response, ...
                             feat_name_resp, ...
                             feat_name_stim, ...
                             harm_str, ...
                             jsonencode(concat_ep), ...
                             snr_str ...
                             ); 

    tbl = readtable(fullfile(par.rsa_path, subdirs, fname_tbl)); 


    load_path = fullfile(par.rsa_path, subdirs); 

    fname = fullfile(load_path, ...
        sprintf('feat-%s_method-pearson_rdm.mat', feat_name_resp)); 

    rdm = load(fname); 
    res.rdm = rdm.rdm; 

    fname = fullfile(load_path, ...
        sprintf('feat-%s_featStim-%s_method-pearson_nCateg-2_perm.mat', ...
                feat_name_resp, feat_name_stim)); 

    perm_res = load(fname); 

    res.best_model = perm_res.best_model; 
    res.boundary = perm_res.all_models_info{1, 'bound1'}; 
    res.r = perm_res.r_obs; 
    res.p_perm = perm_res.p_perm; 

    [res.tested_bounds, idx] = sort(perm_res.all_models_info.bound1); 
    res.tested_bounds = ensure_row(res.tested_bounds); 

    res.r_all_models = ...
                        ensure_row(perm_res.all_models_info.r(idx)); 


    res.r_null_dist = perm_res.r_null_dist; 

    mask = tbl.n_categ_fit == 2; 
    res.z_snr = tbl{mask, 'z_snr'}; 

    mX_data = load(fullfile(load_path, 'mX.mat')); 


    % save data 
    save([fname_save_base, '_data.mat'], 'res'); 

    %% load stimulus 

    load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 

    rdm_stim = rsa_stim(par, 1/par.cycle_dur * harmonics, ...
                           'feat', feat_name_stim); 


    %% individual plots to assemble Figures for the paper

    cmap_ratios = customcolormap([0, 1], {'#e87b0e', '#6a24a3'}, par.n_cond);

    % RHO of individual subjects 
    f = figure('pos', [809 696 80 173], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select(); 
    plot_r_bar([res.r], [res.p_perm], 'y_maxlim', 1, 'ax', ax); 
    pnl.margin = [8,4,0,1]; 

    saveas(f, [fname_save_base, '_rhoInd.svg']); 
    close(f); 

    % overlay of best fitting models across all significant subjects 
    rdms = {res.best_model}; 
    mask_signif = [res.p_perm] < 0.05; 

    f = plot_rdms_overlay(rdms(mask_signif), 'alpha', 1); 

    saveas(f, [fname_save_base, '_bestModelOverlaySignif.svg']); 
    close(f); 



    % grand avearge RDM
    f = plot_rdm(res.rdm); 

    saveas(f, [fname_save_base, '_grandAvgRDM.svg']); 
    close(f); 


    %%

    % individual level permutation test with fixed model
    if exist([fname_save_base, '_permIndFixModel.mat']) && ~run_stats
        clear res_ind_fix_model
        load([fname_save_base, '_permIndFixModel.mat']); 
    else
        res_ind_fix_model = []; 
        for i_sub=1:length(res)
            i_sub
            [res_ind_fix_model(i_sub).r_obs, ...
             res_ind_fix_model(i_sub).r_null_dist, ...
             res_ind_fix_model(i_sub).p_perm, ...
             res_ind_fix_model(i_sub).models, ...
             res_ind_fix_model(i_sub).models_info, ...
             ] = ...
                                test_each_categ_model(...
                                          res(i_sub).rdm, 2, par.n_perm_rsa_group, ...
                                          'covariate', rdm_stim, ...
                                          'n_mini_categ_skip', 1); 
        end
        save([fname_save_base, '_permIndFixModel.mat'], 'res_ind_fix_model'); 
    end


    %%

    close all




end



