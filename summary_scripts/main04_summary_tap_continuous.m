% Summarize RSA results for continuous tapping responses. 

clear 

%% 

experiment = '2ioi'; 

par = get_par('experiment', experiment); 

subjects = num2cell(par.subjects); 


feat_name_resp = 'X'; % X or mX
feat_name_stim = 'ratio';

response = 'tap-force'; % tap-force, tap-impulse

concat_ep = false; 
feat_snr = false; 

harm_selections = {
    [1:12]
}; 

run_stats = false; 

plot_big_figure = true; 

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

    res = []; 

    for i_sub=1:length(subjects)

        res(i_sub).sub = subjects{i_sub}; 
        sub_str = sub_num2str(res(i_sub).sub); 

        load_path_sub = fullfile(par.rsa_path, subdirs, sub_str); 

        fname = fullfile(load_path_sub, ...
            sprintf('%s_feat-%s_method-pearson_rdm.mat', sub_str, feat_name_resp)); 

        rdm = load(fname); 
        res(i_sub).rdm = rdm.rdm; 


        fname = fullfile(load_path_sub, ...
            sprintf('%s_feat-%s_featStim-%s_method-pearson_nCateg-2_perm.mat', ...
                    sub_str, feat_name_resp, feat_name_stim)); 

        perm_res = load(fname); 

        res(i_sub).best_model = perm_res.best_model; 
        res(i_sub).boundary = perm_res.all_models_info{1, 'bound1'}; 
        res(i_sub).r = perm_res.r_obs; 
        res(i_sub).p_perm = perm_res.p_perm; 

        [res(i_sub).tested_bounds, idx] = sort(perm_res.all_models_info.bound1); 
        res(i_sub).tested_bounds = ensure_row(res(i_sub).tested_bounds); 

        res(i_sub).r_all_models = ...
                            ensure_row(perm_res.all_models_info.r(idx)); 


        res(i_sub).r_null_dist = perm_res.r_null_dist; 

        mask = tbl.sub == res(i_sub).sub & tbl.n_categ_fit == 2; 
        res(i_sub).z_snr = tbl{mask, 'z_snr'}; 

    end

    % save data 
    save([fname_save_base, '_data.mat'], 'res'); 

    %% load stimulus 

    load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 

    rdm_stim = rsa_stim(par, 1/par.cycle_dur * harmonics, ...
                           'feat', feat_name_stim); 

    %% save table 

    tbl = table(); 
    tbl.sub = [1 : length(subjects)]'; 
    tbl.zSNR = round([res.z_snr]', 2);
    tbl.boundary = round(bound_to_ratio(stim.ratios(:,1), [res.boundary]), 3);
    tbl.r = round([res.r]', 2);
    tbl.p = [res.p_perm]';

    writetable(tbl, [fname_save_base, '_data.csv']); 

    %% group-level permutation test with free model 

    if exist([fname_save_base, '_permGroup.mat']) && ~run_stats
        clear res_group
        load([fname_save_base, '_permGroup.mat']); 
    else
        res_group = []; 
        [res_group.r_obs, res_group.r_null_dist, res_group.p_perm] = ...
                                        perm_rdm_corr_group(...
                                                  {res.rdm}, 2, 10000, ...
                                                  'covariate', rdm_stim, ...
                                                  'n_mini_categ_skip', 1); 

        save([fname_save_base, '_permGroup.mat'], 'res_group'); 
    end                                        

    %% grand average RDM 

    load_path = fullfile(par.rsa_path, subdirs, 'sub-grand'); 
    fname = sprintf('sub-grand_feat-%s_method-pearson_rdm.mat', feat_name_resp); 
    res_grand = load(fullfile(load_path, fname)); 


    %% grand average mX 

    load_path = fullfile(par.rsa_path, subdirs, 'sub-grand'); 
    fname = 'sub-grand_mX.mat'; 
    mX_data = load(fullfile(load_path, fname)); 



    %% big summary figure 

    if plot_big_figure

        f = plot_big_summary_figure(par, res, res_group, res_grand.rdm, mX_data)

        print(f, '-dpng', '-painters', '-r600', [fname_save_base, '_summary.png']); 

    end




    %% individual plots to assemble Figures for the paper
    cmap_ratios = customcolormap([0, 1], {'#e87b0e', '#6a24a3'}, par.n_cond);

    %% individual RDMs

    f = figure('pos',[809 69 113 959], 'color', 'white'); 
    pnl = panel(f); 
    plot_all_rdms(pnl, ...
        [res.sub], ...
        {res.rdm}, ...
        {res.best_model}, ...
        [res.p_perm] ...
        )
    pnl.margin = [10, 1, 1, 10]; 

    print(f, [fname_save_base, '_indRDMs.eps'], '-depsc', '-painters')
    close(f); 

    %% FIGURE 2

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

    f = plot_rdms_overlay(rdms(mask_signif), 'alpha', 0.1); 

    saveas(f, [fname_save_base, '_bestModelOverlaySignif.svg']); 
    close(f); 


    % median best-fitting model 
    median_bound = median([res.boundary]); 
    median_model = get_2categ_model(13, median_bound); 
    f = figure('pos', [809 756 145 113], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select(); 
    plot_rdm(median_model, 'ax', ax); 
    saveas(f, [fname_save_base, '_medianBestModelRDM.svg']); 
    close(f); 

    rsa_bounds = [res.boundary]; 
    median(rsa_bounds)
    m = bootstrp(1000, @median, rsa_bounds);
    prctile(m, [2.5, 100-2.5])

    bounds_ratios = bound_to_ratio(stim.ratios(:, 1), [res.boundary]); 
    median(bounds_ratios)
    m = bootstrp(1000, @median, bounds_ratios);
    prctile(m, [2.5, 100-2.5])




    % grand avearge RDM
    f = figure('pos', [809 756 145 113], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select(); 
    plot_rdm(res_grand.rdm, 'ax', ax); 
    pnl.margin = [1,1,1,1];

    saveas(f, [fname_save_base, '_grandAvgRDM.svg']); 
    close(f); 


    % MDS            
    cols = customcolormap([0, 1], {'#e87b0e', '#6a24a3'}, par.n_cond);
    Y = {}; 
    for i_sub=1:length(res)
        D = -res(i_sub).rdm + 1; 
        Y{i_sub} = cmdscale(D, 2);
    end
    Y = mean(cat(3, Y{:}), 3); 
    f = figure('pos', [606 688 165 181], 'color', 'white'); 
    s = scatter(Y(:, 1), Y(:, 2), 100, cols, 'filled'); 
    axis([min(Y(:,1)) - 0.2 * range(Y(:,1)), ...
          max(Y(:,1)) + 0.2 * range(Y(:,1)), ...
          min(Y(:,2)) - 0.2 * range(Y(:,2)), ...
          max(Y(:,2)) + 0.2 * range(Y(:,2))])
    axis square           
    ax = gca; 
    ax.XTick = []; 
    ax.YTick = []; 
    ax.XDir = 'reverse'; 

    saveas(f, [fname_save_base, '_MDS.svg']); 
    close(f); 



    %% FIGURE 4 - FIXED MODELS


    % group-level permutation test with fixed model 
    if exist([fname_save_base, '_permGroupFixModel.mat']) && ~run_stats
        clear res_group
        load([fname_save_base, '_permGroupFixModel.mat']); 
    else
        res_group = []; 
        [res_group.r_obs, res_group.r_null_dist, res_group.p_perm, res_group.models, res_group.models_info] = ...
                                perm_rdm_corr_group_fix_model(...
                                          {res.rdm}, 2,  par.n_perm_rsa_group, ...
                                          'covariate', rdm_stim, ...
                                          'n_mini_categ_skip', 1); 

        save([fname_save_base, '_permGroupFixModel.mat'], 'res_group'); 
    end

    % individual level permutation test with fixed model
    if exist([fname_save_base, '_permIndFixModel.mat']) && ~run_stats
        clear res_ind_fix_model
        load([fname_save_base, '_permIndFixModel.mat']); 
    else
        res_ind_fix_model = []; 
        for i_sub=1:length(subjects)
            i_sub
            [res_ind_fix_model(i_sub).r_obs, ...
             res_ind_fix_model(i_sub).r_null_dist, ...
             res_ind_fix_model(i_sub).p_perm, ...
             res_ind_fix_model(i_sub).models, ...
             res_ind_fix_model(i_sub).models_info, ...
             ] = ...
                                test_each_categ_model(...
                                          res(i_sub).rdm, 2,  par.n_perm_rsa_group, ...
                                          'covariate', rdm_stim, ...
                                          'n_mini_categ_skip', 1); 
        end
        save([fname_save_base, '_permIndFixModel.mat'], 'res_ind_fix_model'); 
    end



    f = figure('pos', [652 621 300 311], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select(); 

    n_sub = length(res_ind_fix_model); 
    n_models = length(res_group.r_obs); 

    tmp = {res_ind_fix_model.r_obs}; 
    r_ind = cat(2, tmp{:})';  

    tmp = {res_ind_fix_model.p_perm}; 
    p_ind = cat(1, tmp{:}); 

    c_group = [49, 97, 127] / 255; 
    c_ind = [135, 175, 201]/255; 

    for i_model=1:n_models

        h = bar(ax, i_model, res_group.r_obs(i_model), ...
                'FaceColor', 'none', ...
                'EdgeColor', cmap_ratios(i_model, :), ...
                'LineWidth', 3); 

        hold(ax, 'on'); 

        idx_signif = p_ind(:, i_model) < 0.05; 

        h_signif = scatter(ax, i_model + 0.1 * (rand(1, sum(idx_signif))-0.5), ...
                           r_ind(idx_signif, i_model), ...
                           50, cmap_ratios(i_model, :), 'filled'); 

        h_nonsignif = scatter(ax, i_model + 0.1 * (rand(1, sum(~idx_signif))-0.5), ...
                              r_ind(~idx_signif, i_model), ...
                              50, cmap_ratios(i_model, :)); 

        h_signif.MarkerFaceAlpha = 0.6;                   
        h_nonsignif.MarkerFaceAlpha = 0.6;                   

        r_sem = std(r_ind(:, i_model)) / sqrt(n_sub); 
        r_ci = norminv(1 - 0.025) * r_sem; 

        plot([i_model, i_model], ...
             [res_group.r_obs(i_model) - r_ci, ...
              res_group.r_obs(i_model) + r_ci], ...
              'color', 'k', 'linew', 4)

        txt = text(i_model, max(r_ind(:,i_model))+0.1, sprintf('%.2g', res_group.p_perm(i_model))); 
        txt.HorizontalAlignment = 'center'; 

    end




    ax.YLim = [-1, 1]; 
    ax.YTick = [-1, 0, 1]; 
    ax.YAxis.TickDirection = 'out'; 
    ax.XAxis.TickDirection = 'out'; 

    ax.XTick = [1 : n_models]; 
    ax.XTickLabel = sort(res_group.models_info.bound1); 


    saveas(f, [fname_save_base, '_fixedModels.svg']); 

    %%

    f = plot_r_over_all_models(cat(1, res.r_all_models), res(1).tested_bounds);

    saveas(f, [fname_save_base, '_fixedModelsDistr.svg']); 


    %%

    f = figure('pos',[652 45 175 887], 'color', 'white'); 
    pnl = panel(f); 
    pnl.pack('v', length(res)); 
    for i_sub=1:length(res)
        ax = pnl(i_sub).select(); 
        plot_r_over_all_models(cat(1, res(i_sub).r_all_models), res(1).tested_bounds, 'ax', ax);
        ax.XTickLabel = []; 
        ax.XAxis.Label.String = ''; 
        ax.YAxis.Label.String = ''; 
    end
    pnl.de.margin = [1, 1, 1, 0]; 

    saveas(f, [fname_save_base, '_fixedModelsDistrInd.svg']); 


    %%

    close all




end

