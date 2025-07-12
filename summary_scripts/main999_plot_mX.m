% Plots grand average magnitde spectra for EEG, tap-force, and tap-impulse
% signals. 

clear 

%% 

experiment = '2ioi'; 

subjects = num2cell([10, 11]); 

response = 'eeg'; 

roi_name = 'all';         

concat_ep = false; 
feat_snr = false; 


%%

par = get_par('experiment', experiment); 


if feat_snr 
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 
else
    snr_str = 'none'; 
end

save_path = fullfile(par.deriv_path, 'figures', 'fig_mX'); 
mkdir(save_path); 

fname_save_base = fullfile(save_path, ...
    sprintf('roi-%s_snr-%s', ...
                     roi_name, ...
                     snr_str ...
                     ));


%% STIM 

load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 

env = cellfun(@(x) abs(hilbert(x)), stim.s, 'uni', 0); 

[X, freq_stim] = get_X(env, stim.fs); 
mX_stim = cellfun(@(x) abs(x), X, 'uni', 0); 

harmonics = [1 : floor(freq_stim(end) / (1/par.cycle_dur))]; 
frex_stim = 1/par.cycle_dur * harmonics; 

%% TAPPING

mX_all = cell(1, length(subjects)); 
for i_sub=1:length(subjects)
    sub = subjects{i_sub}; 
    sub_str = sub_num2str(sub);   
    % load data 
    fpath = fullfile(par.deriv_path, 'fft', sub_str); 
    fname = sprintf('sub-%03d_response-tap-impulse_concatEp-false_snr-none_mX', sub); 
    [header, data] = CLW_load(fullfile(fpath, fname));
    mX_all{i_sub} = data; 
end

% average mX across all subjects 
mX_tap_impulse = squeeze(mean(cat(7, mX_all{:}), 7)); 
freq_tap_impulse = [0 : header.datasize(end)-1] * header.xstep; 

%     harmonics = [1 : floor(freq_tap_impulse(end) / (1/par.cycle_dur))]; 
harmonics = [1 : 12]; 
frex_tap_impulse = 1/par.cycle_dur * harmonics; 

%% TAPPING

mX_all = cell(1, length(subjects)); 
for i_sub=1:length(subjects)
    sub = subjects{i_sub}; 
    sub_str = sub_num2str(sub);   
    % load data 
    fpath = fullfile(par.deriv_path, 'fft', sub_str); 
    fname = sprintf('sub-%03d_response-tap-force_concatEp-false_snr-none_mX', sub); 
    [header, data] = CLW_load(fullfile(fpath, fname));
    mX_all{i_sub} = data; 
end

% average mX across all subjects 
mX_tap_force = squeeze(mean(cat(7, mX_all{:}), 7)); 
freq_tap_force = [0 : header.datasize(end)-1] * header.xstep; 

%     harmonics = [1 : floor(freq_tap_force(end) / (1/par.cycle_dur))]; 
harmonics = [1 : 12]; 
frex_tap_force = 1/par.cycle_dur * harmonics; 

%% EEG

mX_all = cell(1, length(subjects)); 
for i_sub=1:length(subjects)
    sub = subjects{i_sub}; 
    sub_str = sub_num2str(sub);   
    % load data 
    fpath = fullfile(par.deriv_path, 'fft', sub_str); 
    fname = sprintf('sub-%03d_response-eeg_concatEp-false_snr-2-5_mX', sub); 
    [header, data] = CLW_load(fullfile(fpath, fname));
    % average frontocentral 
    if strcmp(roi_name, 'all')
        chan_idx = [1 : length(header.chanlocs)];  
    elseif strcmp(roi_name, 'front')            
        chan_idx = get_chan_idx(header, par.roi_front_chans); 
    end
    data = mean(data(:, chan_idx, 1, 1, 1, :), 2); 
    mX_all{i_sub} = abs(data); 
end

% average mX across all subjects 
mX_eeg = squeeze(mean(cat(7, mX_all{:}), 7)); 
freq_eeg = [0 : header.datasize(end)-1] * header.xstep; 

%     harmonics = [1 : floor(freq_eeg(end) / (1/par.cycle_dur))]; 
harmonics = [1 : 12]; 
frex_eeg = 1/par.cycle_dur * harmonics; 

%%

f = figure('pos', [652 150 740 782], 'color', 'white'); 
pnl = panel(f); 
pnl.pack('h', 4); 

plot_mX(mX_stim, freq_stim, frex_stim, 'pnl', pnl(1), 'fmin', 1, 'fmax', 17, 'show_x_line', true); 
plot_mX(num2cell(mX_tap_impulse, 2), freq_tap_impulse, frex_tap_impulse, 'pnl', pnl(2), 'fmin', 1, 'fmax', 17, 'show_x_line', true); 
plot_mX(num2cell(mX_tap_force, 2), freq_tap_force, frex_tap_force, 'pnl', pnl(3), 'fmin', 1, 'fmax', 17, 'show_x_line', true); 
plot_mX(num2cell(mX_eeg, 2), freq_eeg, frex_eeg, 'pnl', pnl(4), 'fmin', 1, 'fmax', 17, 'show_x_line', true); 

print(f, [fname_save_base, '_mX.eps'], '-depsc', '-painters')


