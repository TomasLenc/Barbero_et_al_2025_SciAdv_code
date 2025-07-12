% Use tap onset times to calculate produced inter-tap interval ratios.  

clear 

par = get_par('experiment', '2ioi'); 

sub = 15; 

%%

load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 
stim_ratios = stim.ratios; 
stim_iois = stim.ratios * stim.cycle_dur; 

%% 

sub_str = sub_num2str(sub); 

save_path_onsets = fullfile(par.deriv_path, 'preproc', sub_str); 
mkdir(save_path_onsets); 

par_onset_extract = []; 
par_onset_extract.tap_onset_thr = 0.1; 
par_onset_extract.tap_onset_min_asy = 0.030; 
par_onset_extract.tap_tone_min_asy = 0.080; 
par_onset_extract.n_cycles_skip = 0; 

[header, data] = CLW_load(fullfile(par.raw_path, sprintf('%s_task-TAP.lw6', sub_str))); 

fs = 1/header.xstep; 

tap_onset_times_all = cell(stim.n_cond, 1);
asy_all = cell(stim.n_cond, 1); 

c = 1; 
for i_cond=1:13        
    
    cond_idx = find(strcmp({header.events.code}, num2str(i_cond))); 

    tap_onset_times_cond = cell(1, length(cond_idx)); 
    asy_cond = cell(1, length(cond_idx)); 

    for i_ep=1:length(cond_idx)

        fz = ensure_row(squeeze(data(cond_idx(i_ep), 1, 1, 1, 1, :))); 
        
        s = ensure_row(squeeze(data(cond_idx(i_ep), 2, 1, 1, 1, :))); 
        
        trigger_on = ensure_row(squeeze(data(cond_idx(i_ep), 3, 1, 1, 1, :))); 

        % get tap onset times
        [tap_onset_times, tap_indices] = extract_taps(...
                                       trigger_on, ...
                                       fs, ...
                                       par_onset_extract.tap_onset_min_asy, ...
                                       par_onset_extract.tap_onset_min_asy, ...
                                       'plot_diagnostic', false); 

        % estimtae asynchronies
        asy = get_asynchronies(...
                    tap_onset_times, ...
                    stim_iois(i_cond, :), ...
                    par.trial_dur / par.cycle_dur, ...
                    'n_cycles_skip', par_onset_extract.n_cycles_skip); 

        asy_cond{i_ep} = asy; 
        tap_onset_times_cond{i_ep} = tap_onset_times; 

        c = c+1; 
        close all

    end     

    tap_onset_times_all{i_cond} = tap_onset_times_cond; 
    asy_all{i_cond} = asy_cond; 

end

% get mean asynchrony 
asy_all = cat(2, asy_all{:}); 
asy_all = cat(2, asy_all{:}); 
asy_all(abs(asy_all) > 0.300) = [];     
mean_asy = mean(asy_all); 

tap_ratios_all = cell(stim.n_cond, 1); 
n_valid_tap_ratios_all = zeros(stim.n_cond, 1); 

for i_cond=1:13        

    for i_ep=1:3

        % get ratios 
        [tap_ratios, ~, n_valid_tap_ratios] = ...
                    extract_tap_ratios(...
                            tap_onset_times_all{i_cond}{i_ep}, ...
                            stim_iois(i_cond, :), ...
                            par.trial_dur / par.cycle_dur, ...
                            'tap_tone_min_asy', par_onset_extract.tap_tone_min_asy, ...
                            'mean_asy', mean_asy, ...
                            'n_cycles_skip', par_onset_extract.n_cycles_skip, ...
                            'plot_diagnostic', false); 
                        
            t = [0 : length(trigger_on)-1] / fs; 

%             % diagnostic figure (sanity check)
%             figure('Position', [620 727 1162 320])
%             plot(t, trigger_on, 'linew', 2); 
%             hold on 
%             plot(t, s); 
%             fz_norm = fz; 
%             fz_norm = fz_norm ./ max(fz_norm); 
%             fz_norm = fz_norm - min(fz_norm); 
%             plot(t, fz_norm, 'color', [.5 .5 .5]); 
%             plot([tap_onset_times_all{i_cond}{i_ep}', tap_onset_times_all{i_cond}{i_ep}'], [-1, 2], 'm--')
%             xlim([-1, 23]); 
%             pause

        tap_ratios_all{i_cond} = [tap_ratios_all{i_cond}; tap_ratios]; 
        n_valid_tap_ratios_all(i_cond) = n_valid_tap_ratios_all(i_cond) + n_valid_tap_ratios; 

        close all

    end

end

% get mean tap ratio per condition
tap_ratios_mean = cellfun(@(x) mean(x, 1), tap_ratios_all, 'uni', 0); 
tap_ratios_mean = cat(1, tap_ratios_mean{:}); 

% save tap ratios 
fname = sprintf('sub-%03d_response-tap-ratios.mat', sub); 
tap_ratios = tap_ratios_mean; 
n_valid_tap_ratios = n_valid_tap_ratios_all; 

save(fullfile(save_path_onsets, fname), ...
    'tap_onset_times_all', 'tap_ratios', ...
    'stim_ratios', 'tap_ratios_all', 'n_valid_tap_ratios'); 

