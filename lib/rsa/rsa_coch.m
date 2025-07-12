function rsa_coch(par_rsa, varargin)

parser = inputParser(); 

addParameter(parser, 'response', 'urear-an');
addParameter(parser, 'feat_resp', 'X'); 
addParameter(parser, 'feat_stim', 'X'); 
addParameter(parser, 'rerun_stats', true); 
addParameter(parser, 'rerun_figs', true); 
addParameter(parser, 'par', []); 
addParameter(parser, 'flip_rdm', true); 

parse(parser, varargin{:}); 

response = parser.Results.response; 
feat_name_resp = parser.Results.feat_resp; 
feat_name_stim = parser.Results.feat_stim; 
rerun_stats = parser.Results.rerun_stats; 
rerun_figs = parser.Results.rerun_figs; 
par = parser.Results.par; 
flip_rdm = parser.Results.flip_rdm; 

% if parameter structure hasn't been explicitly passed, load a new one from
% the file - this allows the possibility for the user to pass an adjusted
% par structure intead of the default one
if isempty(par)
    par = get_par(); 
end


%% 
% frequencies of interest (par_rsa.harmonics of rhythm cycle repetition rate) 

harm_str = vec_to_str(par_rsa.harmonics, 'format','%d', 'sep','-'); 

frex = 1 / par.cycle_dur * par_rsa.harmonics; 

%% input directories 

load_path = fullfile(par.deriv_path, 'fft', sprintf('response-%s', response)); 


%% output directories

if par_rsa.feat_snr 
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 
else
    snr_str = 'none'; 
end

data_path = fullfile(par.rsa_path, ...
                         sprintf('response-%s', response), ...
                         sprintf('harm-%s', harm_str), ...
                         sprintf('concatEp-%s_snr-%s', jsonencode(par_rsa.concat_ep), snr_str) ...
                         ); 

save_path = fullfile(data_path); 

mkdir(save_path); 

%% open csv file 

fname_tbl = sprintf('response-%s_feat-%s_featStim-%s_harm-%s_concatEp-%s_snr-%s.csv', ...
                         response, ...
                         feat_name_resp, ...
                         feat_name_stim, ...
                         harm_str, ...
                         jsonencode(par_rsa.concat_ep), ...
                         snr_str); 

if exist(fullfile(data_path, fname_tbl))
    tbl_res = readtable(fullfile(data_path, fname_tbl)); 
else
    col_names = {'n_categ_fit', 'bound1', 'bound2', 'r', 'p', 'p_adj', ...
                'z_snr', 'z_snr_sd', 'n_models_fit', 'n_perm'}; 
    tbl_res = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 
end


%% load stimulus info 

load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 

rdm_stim = rsa_stim(par, frex, 'feat', feat_name_stim); 

%% load FFTs 

% load complex FFT without SNR subtraction 
fname = sprintf('response-%s_concatEp-%s_snr-none_X', ...
            response, jsonencode(par_rsa.concat_ep)); 

[header_X, data_X] = CLW_load(fullfile(load_path, fname)); 

% load complex FFT with SNR subtraction 
fname = sprintf('response-%s_concatEp-%s_snr-%d-%d_X', ...
            response, ...
            jsonencode(par_rsa.concat_ep), ...
            par.snr_bins(1), par.snr_bins(2)); 

[header_X_snr, data_X_snr] = CLW_load(fullfile(load_path, fname)); 

% load magniude FFT without SNR subtraction 
fname = sprintf('response-%s_concatEp-%s_snr-none_mX', ...
            response, jsonencode(par_rsa.concat_ep)); 

[header_mX, data_mX] = CLW_load(fullfile(load_path, fname)); 

% load magniude FFT with SNR subtraction 
fname = sprintf('response-%s_concatEp-%s_snr-%d-%d_mX', ...
            response, jsonencode(par_rsa.concat_ep), ...
            par.snr_bins(1), par.snr_bins(2)); 

[header_mX_snr, data_mX_snr] = CLW_load(fullfile(load_path, fname)); 


%% get z_snr

freq = [0 : header_mX.datasize(end)-1] * header_mX.xstep;   

% get zSNR from raw magnitude spectra
z_snr = get_z_snr(data_mX, freq, frex, par.snr_bins(1), par.snr_bins(2)); 
amp_sum = get_amp_summary(data_mX_snr, freq, frex, 'method', 'sum'); 


%% get RDM 
    
% get features and RDM from complex-valued spectra 
if par_rsa.feat_snr
    [feat_for_rdm] = get_feat_from_X(squeeze(data_X_snr), freq, frex, ...
                                     'feat', feat_name_resp); 
else
    [feat_for_rdm] = get_feat_from_X(squeeze(data_X), freq, frex, ...
                                     'feat', feat_name_resp);    
end

rdm = get_rdm_from_feat(feat_for_rdm, 'method', 'pearson'); 


%%

% % plot magmitude spectra 
fname = sprintf('mX'); 

if par_rsa.feat_snr
    header_mX_to_plot = header_mX_snr; 
    data_mX_to_plot = data_mX_snr; 
else
    header_mX_to_plot = header_mX; 
    data_mX_to_plot = data_mX; 
end

data_to_save = [];
data_to_save.header_mX_to_plot = header_mX_to_plot; 
data_to_save.data_mX_to_plot = data_mX_to_plot; 
data_to_save.frex = frex; 
data_to_save.z_snr = z_snr; 
save(fullfile(save_path, [fname, '.mat']), '-struct', 'data_to_save'); 
 
% if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs
%     
%     f = plot_mX_lw(header_mX_to_plot, data_mX_to_plot, 1, frex, ...
%                 'z_snr', z_snr, 'fmax', max([frex, 15])); 
% 
%     print(f, '-dpng', '-painters', '-r600', fullfile(save_path, fname))
% 
% end

% RDM 
fname = sprintf('feat-%s_method-pearson_rdm', feat_name_resp); 
save(fullfile(save_path, [fname, '.mat']), 'rdm'); 

% if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs
% 
%     f = plot_rdm(rdm, 'ratios', stim.ratios, 'flip', flip_rdm); 
% 
%     print(f, '-dpng', '-painters', '-r300', fullfile(save_path, fname))
% 
%     close all
% end

%% 

n_categ = 2; 

fname = sprintf('feat-%s_featStim-%s_method-pearson_nCateg-2', ...
    feat_name_resp, feat_name_stim); 

fname_mat = fullfile(save_path, [fname, '_perm.mat']); 

if isfile(fname_mat) & ~rerun_stats
    
    warning('loading permutation data from file `%s`', fname); 
    
    perm_output = load(fname_mat); 
    
    r_obs = perm_output.r_obs; 
    p_perm = perm_output.p_perm; 
    r_null_dist = perm_output.r_null_dist; 
    best_model = perm_output.best_model; 
    all_models_info = perm_output.all_models_info; 
    
else
    
    [r_obs, p_perm, r_null_dist, best_model, best_model_info, all_models_info] = ...
                          find_test_best_categ_model(...
                                        rdm, n_categ, par.n_perm_rsa_individual, ...
                                        'covariate', rdm_stim, ...
                                        'n_mini_categ_skip', par.n_mini_categ_skip, ...
                                        'verbose', false);
                                    
    save(fullfile(save_path, [fname, '_perm.mat']), ...
                  'par', 'par_rsa', 'r_obs', 'p_perm', 'r_null_dist', ...
                  'best_model', 'all_models_info');  
                   
    writetable(all_models_info, fullfile(save_path, ...
               [fname, '_models.csv'])); 
end

% write table 
row_mask = (tbl_res.n_categ_fit == 2); 
tbl_res(row_mask, :) = []; 

new_row = [...
    {2}, ...
    all_models_info{1, 'bound1'}, ...
    all_models_info{1, 'bound2'}, ...
    all_models_info{1, 'r'}, ...
    all_models_info{1, 'p'}, ...
    all_models_info{1, 'p_adj'}, ...
    mean(z_snr), ...
    std(z_snr), ...
    size(all_models_info, 1), ...
    par.n_perm_rsa_individual, ...
    ]; 

tbl_res = [tbl_res; new_row]; 

% if ~isfile(fullfile(save_path, [fname, '_perm.png'])) || rerun_figs
% 
%     f = plot_best_model(rdm, best_model, r_obs, p_perm, r_null_dist, ...
%                         'z_snr', z_snr, 'ratios', stim.ratios, 'flip', flip_rdm);
% 
%     print(f, '-dpng', '-painters', '-r300', ...
%           fullfile(save_path, [fname, '_perm']))
% 
%     close all
% end


%% 

writetable(tbl_res, fullfile(data_path, fname_tbl)); 

