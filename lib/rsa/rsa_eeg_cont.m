function rsa_eeg_cont(sub, par_rsa, varargin)

parser = inputParser(); 

addParameter(parser, 'feat_resp', 'X'); 
addParameter(parser, 'feat_stim', 'ratio'); % X, ratio, urear-an
addParameter(parser, 'rerun_stats', true); 
addParameter(parser, 'rerun_figs', true); 
addParameter(parser, 'par', []); 
addParameter(parser, 'flip_rdm', true); 

parse(parser, varargin{:}); 

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

sub_str = sub_num2str(sub); 

load_path = fullfile(par.deriv_path, 'fft', sub_str); 


%% output directories

if par_rsa.feat_snr 
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 
else
    snr_str = 'none'; 
end

subdirs = fullfile('response-eeg', ...
                     sprintf('harm-%s', harm_str), ...
                     sprintf('concatEp-%s_snr-%s', jsonencode(par_rsa.concat_ep), snr_str), ...
                     sprintf('roi-%s_avgChans-%s', par_rsa.roi_name, jsonencode(par_rsa.avg_chans)) ...
                     ); 

save_path = fullfile(par.rsa_path, subdirs, sub_str); 

mkdir(save_path); 
    
%% open csv file 

fname_tbl = sprintf('response-eeg_feat-%s_featStim-%s_harm-%s_concatEp-%s_snr-%s_roi-%s_avgChans-%s.csv', ...
                         feat_name_resp, ...
                         feat_name_stim, ...
                         harm_str, ...
                         jsonencode(par_rsa.concat_ep), ...
                         snr_str, ...
                         par_rsa.roi_name, ...
                         jsonencode(par_rsa.avg_chans) ...
                         ); 

if exist(fullfile(par.rsa_path, subdirs, fname_tbl))
    tbl_res = readtable(fullfile(par.rsa_path, subdirs, fname_tbl)); 
else
    col_names = {'sub', 'n_categ_fit', 'bound1', 'bound2', 'r', 'p', 'p_adj', ...
                'z_snr', 'z_snr_sd', 'n_models_fit', 'n_perm', 'chans'}; 
    tbl_res = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 
end


%% load stimulus info 

if strcmp(feat_name_stim, 'urear-an')
    coch_path = fullfile(par.rsa_path, ...
                             sprintf('response-%s', feat_name_stim), ...
                             sprintf('harm-%s', harm_str), ...
                             sprintf('concatEp-%s_snr-%s', jsonencode(par_rsa.concat_ep), snr_str) ...
                             ); 
    fname = sprintf('feat-%s_method-pearson_rdm', feat_name_resp); 
    coch = load(fullfile(coch_path, [fname, '.mat']), 'rdm');    
    rdm_stim = coch.rdm; 
else
    load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 
    rdm_stim = rsa_stim(par, frex, 'feat', feat_name_stim); 
end

%% load FFTs 

% load complex FFT without SNR subtraction 
fname = sprintf('%s_response-eeg_concatEp-%s_snr-none_X', ...
            sub_str, jsonencode(par_rsa.concat_ep)); 

[header_X, data_X] = CLW_load(fullfile(load_path, fname)); 

% load complex FFT with SNR subtraction 
fname = sprintf('%s_response-eeg_concatEp-%s_snr-%d-%d_X', ...
            sub_str, jsonencode(par_rsa.concat_ep), ...
            par.snr_bins(1), par.snr_bins(2)); 

[header_X_snr, data_X_snr] = CLW_load(fullfile(load_path, fname)); 

% load magniude FFT without SNR subtraction 
fname = sprintf('%s_response-eeg_concatEp-%s_snr-none_mX', ...
            sub_str, jsonencode(par_rsa.concat_ep)); 

[header_mX, data_mX] = CLW_load(fullfile(load_path, fname)); 

% load magniude FFT with SNR subtraction 
fname = sprintf('%s_response-eeg_concatEp-%s_snr-%d-%d_mX', ...
            sub_str, jsonencode(par_rsa.concat_ep), ...
            par.snr_bins(1), par.snr_bins(2)); 

[header_mX_snr, data_mX_snr] = CLW_load(fullfile(load_path, fname)); 

% save a copy of all channels chanlocs for topoplot 
chanlocs_all_chans = header_X.chanlocs; 

%% get z_snr

% We first need to get the measure of SNR separately for each channel. This
% is needed for (1) topoplots, (2) to sometimes select the most responsoive
% channels for feature extraction. 

freq = [0 : header_mX.datasize(end)-1] * header_mX.xstep;   

z_snr_all_chan = get_z_snr(data_mX, freq, frex, par.snr_bins(1), par.snr_bins(2)); 

amp_sum_all_chan = get_amp_summary(data_mX_snr, freq, frex, 'method', 'sum'); 



%% channel selection

if contains(par_rsa.roi_name, 'best')

    % average across all conditions, separately for each channel (this will be
    % needed to compute best channel)
    amp_sum_mean_cond = mean(amp_sum_all_chan, 1);
    
    % how many most-responsive chanels to use? 
    tmp = regexp(par_rsa.roi_name, '^\d+', 'match'); 
    n_best_chans_to_use = str2num(tmp{:}); 

    % sort the channels 
    [~, best_chan_idx] = sort(amp_sum_mean_cond, 'descend');

    % pick N most overall responsive channels 
    roi_chans_idx = best_chan_idx(1:n_best_chans_to_use); 

    roi_chans = {header_X.chanlocs(roi_chans_idx).labels}; 
    
    fname = sprintf('%s_feat-ampSum_bestChans.csv', sub_str); 
    
    rows = [ensure_col({header_X.chanlocs(best_chan_idx).labels}), ...
            ensure_col(num2cell(amp_sum_mean_cond(best_chan_idx)))]; 

    tbl_best_chan = cell2table(rows, 'VariableNames', {'label', 'ampSum'}); 
    
    writetable(tbl_best_chan, fullfile(save_path, fname)); 

    roi_name_for_tbl = join(roi_chans, '-'); 
    
elseif strcmp(par_rsa.roi_name, 'front')

    % ROI for audio
    roi_chans = {
        'F1'
        'Fz'
        'F2'
        'FC1'
        'FCz'
        'FC2'
        'C1'
        'Cz'
        'C2'
    }; 

    roi_chans_idx = get_chan_idx(header_X, roi_chans); 

    roi_name_for_tbl = {'front'}; 
                           
elseif strcmp(par_rsa.roi_name, 'all')

    % ROI for audio
    roi_chans = {header_X.chanlocs.labels};

    % remove masoids if present 
    mask_mast = ismember(roi_chans, {'TP9', 'TP10'}); 
    if sum(mask_mast) > 2
        error('found more than 2 mastoids in the data??? weird...stopping...')
    end
    roi_chans(mask_mast) = []; 
    
    roi_chans_idx = [1 : length(roi_chans)];
                           
    roi_name_for_tbl = {'all'}; 
    
elseif all(ismember(strsplit(par_rsa.roi_name, '-'), {header_X.chanlocs.labels}))

    roi_chans = strsplit(par_rsa.roi_name, '-'); 
    
    roi_chans_idx = get_chan_idx(header_X, roi_chans); 

    roi_name_for_tbl = join(roi_chans, '-'); 
    
else
   
    error('roi name %s not recognized', par_rsa.roi_name); 
    
end

roi_chans


%% select (and average) channels of interest 

[header_X, data_X] = RLW_arrange_channels(header_X, data_X, roi_chans); 
[header_mX, data_mX] = RLW_arrange_channels(header_mX, data_mX, roi_chans); 
[header_X_snr, data_X_snr] = RLW_arrange_channels(header_X_snr, data_X_snr, roi_chans); 
[header_mX_snr, data_mX_snr] = RLW_arrange_channels(header_mX_snr, data_mX_snr, roi_chans); 

if par_rsa.avg_chans
    
    [header_X, data_X] = RLW_pool_channels(header_X, data_X, roi_chans, ...
                                           'mixed_channel_label', roi_name_for_tbl{:}, ...
                                           'keep_original_channels', false); 
    
    [header_X_snr, data_X_snr] = RLW_pool_channels(header_X_snr, data_X_snr, roi_chans, ...
                                           'mixed_channel_label', roi_name_for_tbl{:}, ...
                                           'keep_original_channels', false); 
    
    % if the parameter is set to averaging channels, it means that the
    % complex DFT across the ROI channels will be averaged before features
    % are extracted for RDM. Now, that is equivalent to avergaging
    % time-domain across those channels, and so there may be some
    % cancellation due to phase differences. In order to make the z-snr and
    % magnitude spectrum plots as useful as possible, let's re-estimate the
    % magnitude spectra from the channel-aveaged complex spectra. 
    header_mX = header_X; 
    data_mX = abs(data_X); 
    
    header_mX_snr = header_X_snr; 
    data_mX_snr = abs(data_X_snr); 
        
end

% re-estimate z-SNR 
z_snr_roi = get_z_snr(data_mX, freq, frex, par.snr_bins(1), par.snr_bins(2)); 
z_snr_roi = mean(z_snr_roi, 2); 


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


%% plotting 

% plot topo 
% ---------

% if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs
% 
%     % prepare custom colormap for topoplot (interpolates between 2 colors - nice!)
%     % cmap = brewermap(64, 'Reds'); 
%     % cmap = customcolormap([0, 1], { '#eb4034', '#ffffff'}, 64);
% 
%     cmap = parula(); 
% 
%     % amp sum
%     % -------
% 
%     fname = sprintf('%s_feat-ampSum_topo', sub_str); 
%     
%     % plot overall response magnitude disstribution across scalp 
%     f = plot_topoplots(amp_sum_all_chan, chanlocs_all_chans, ...
%                   'lab', ...
%                   'amp-sum', ...
%                   'mark_chan_idx', roi_chans_idx, ...
%                   'cmap', cmap); 
% 
% 
%     print(f, '-dpng', '-painters', '-r90', fullfile(save_path, fname));
% 
%     % z SNR
%     % -----
% 
%     fname = sprintf('%s_feat-zSNR_topo', sub_str); 
% 
%     f = plot_topoplots(z_snr_all_chan, chanlocs_all_chans, ...
%                   'lab', ...
%                   'z-snr', ...
%                   'mark_chan_idx', roi_chans_idx, ...
%                   'cmap', cmap); 
% 
% 
%     print(f, '-dpng', '-painters', '-r90', fullfile(save_path, fname));
% 
%     close all
%     
%     
%     close all
% end


% plot magmitude spectra 
% ----------------------
fname = sprintf('%s_mX', sub_str); 

if par_rsa.feat_snr
    header_mX_to_plot = header_mX_snr; 
    data_mX_to_plot = data_mX_snr; 
else
    header_mX_to_plot = header_mX; 
    data_mX_to_plot = data_mX; 
end

if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs

    f = plot_mX_lw(header_mX_to_plot, data_mX_to_plot, ...
                   [1 : header_mX_to_plot.datasize(2)], ...
                   frex, ...
                   'z_snr', z_snr_roi, ...
                   'fmax', max([frex, 15])); 

    print(f, '-dpng', '-painters', '-r200', fullfile(save_path, fname))

end

data_to_save = []; 
data_to_save.header_mX_to_plot = header_mX_to_plot; 
data_to_save.data_mX_to_plot = data_mX_to_plot; 
data_to_save.frex = frex; 
data_to_save.z_snr = z_snr_roi; 

save(fullfile(save_path, [fname, '.mat']), '-struct', 'data_to_save'); 


% plot RDM 
% --------
fname = sprintf('%s_feat-%s_method-pearson_rdm', sub_str, feat_name_resp); 

if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs

    f = plot_rdm(rdm, 'ratios', stim.ratios, 'flip', flip_rdm); 

    print(f, '-dpng', '-painters', '-r200', fullfile(save_path, fname))

    close all

end

save(fullfile(save_path, [fname, '.mat']), 'rdm'); 


%% 2 category models

n_categ = 2; 

fname = sprintf('%s_feat-%s_featStim-%s_method-pearson_nCateg-2', sub_str, ...
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

% write to table 
if isnumeric(sub)
    row_mask = (tbl_res.sub == sub) & (tbl_res.n_categ_fit == 2); 
    tbl_res(row_mask, :) = []; 

    new_row = [...
        {sub}, ...
        2, ...
        all_models_info{1, 'bound1'}, ...
        all_models_info{1, 'bound2'}, ...
        all_models_info{1, 'r'}, ...
        all_models_info{1, 'p'}, ...
        all_models_info{1, 'p_adj'}, ...
        mean(z_snr_roi), ...
        std(z_snr_roi), ...
        size(all_models_info, 1), ...
        par.n_perm_rsa_individual, ...
        roi_name_for_tbl, ...
        ]; 

    tbl_res = [tbl_res; new_row]; 
end


%% 

writetable(tbl_res, fullfile(par.rsa_path, subdirs, fname_tbl)); 

