% This script performs trial averaging and FFT. 

clear

par = get_par('experiment', '2ioi'); 

sub = 15; 

%% 

% prepare output directories
sub_str = sub_num2str(sub); 

save_path = fullfile(par.deriv_path, 'fft', sub_str); 

mkdir(save_path); 


%% low-frequency response 

% load data 
fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt*.lw6', sub); 
d = dir(fullfile(fpath, fname)); 

[header_lf_ep, data_lf_ep] = CLW_load(fullfile(d.folder, d.name));

% re-segment 
codes = unique({header_lf_ep.events.code}); 

[header_lf_ep, data_lf_ep] = segment_safe(...
                                        header_lf_ep, ...
                                        data_lf_ep, ...
                                        codes,...
                                        'x_start', 0, ...
                                        'x_duration', par.trial_dur); 
                                    
% downsample  
[header_lf_ep, data_lf_ep] = RLW_downsample(header_lf_ep, data_lf_ep, ...
                                'x_downsample_ratio', par.decim_factor); 

% rereference 
[header_lf_ep, data_lf_ep] = RLW_rereference(header_lf_ep, data_lf_ep, ...
                                'apply_list', {header_lf_ep.chanlocs.labels}, ...
                                'reference_list', par.eeg_ref_chans); 
    

% don't concatenate epochs
[header_avg, data_avg] = average_trials_per_cond(header_lf_ep, data_lf_ep, ...
                            'concatenate_epochs', false); 

fname = sprintf('sub-%03d_response-eeg_concatEp-false', sub); 

get_and_save_ffts(header_avg, data_avg, par.snr_bins, ...
                  'fname', fname, 'save_path', save_path); 

% concatenate epochs
[header_avg, data_avg] = average_trials_per_cond(header_lf_ep, data_lf_ep, ...
                            'concatenate_epochs', true); 

fname = sprintf('sub-%03d_response-eeg_concatEp-true', sub); 

get_and_save_ffts(header_avg, data_avg, par.snr_bins, ...
                  'fname', fname, 'save_path', save_path); 
                            
                            

%% tapping response 

tap_response_type = 'force'; 

% load data 
fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
fname =  sprintf('sub-%03d_task-TAP_response-%s.lw6', sub, tap_response_type); 
d = dir(fullfile(fpath, fname)); 

[header_tap_env, data_tap_env] = CLW_load(fullfile(d.folder, d.name));

% segment all trials
codes = unique({header_tap_env.events.code}); 

[header_tap_env_ep, data_tap_env_ep] = segment_safe(...
                                        header_tap_env, ...
                                        data_tap_env, ...
                                        codes,...
                                        'x_start', 0, 'x_duration', ...
                                        par.trial_dur); 
                                    
% don't concatenate epochs
[header_avg, data_avg] = average_trials_per_cond(...
                            header_tap_env_ep, data_tap_env_ep, ...
                            'concatenate_epochs', false); 

fname = sprintf('sub-%03d_response-tap-force_concatEp-false', sub); 

get_and_save_ffts(header_avg, data_avg, par.snr_bins, ...
                  'fname', fname, 'save_path', save_path, 'fmax', 100); 

% concatenate epochs
[header_avg, data_avg] = average_trials_per_cond(...
                            header_tap_env_ep, data_tap_env_ep, ...
                            'concatenate_epochs', true); 

fname = sprintf('sub-%03d_response-tap-force_concatEp-true', sub); 

get_and_save_ffts(header_avg, data_avg, par.snr_bins, ...
                  'fname', fname, 'save_path', save_path, 'fmax', 100); 


%% tapping impulses  

tap_response_type = 'impulse'; 

% load data 
fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
fname =  sprintf('sub-%03d_task-TAP_response-%s.lw6', sub, tap_response_type); 
d = dir(fullfile(fpath, fname)); 

[header_tap_imp, data_tap_imp] = CLW_load(fullfile(d.folder, d.name));

% segment all trials
codes = unique({header_tap_imp.events.code}); 

[header_tap_imp, data_tap_imp] = segment_safe(...
                                        header_tap_imp, ...
                                        data_tap_imp, ...
                                        codes,...
                                        'x_start', 0, 'x_duration', ...
                                        par.trial_dur); 
                                    
% don't concatenate epochs
[header_avg, data_avg] = average_trials_per_cond(...
                            header_tap_imp, data_tap_imp, ...
                            'concatenate_epochs', false); 

fname = sprintf('sub-%03d_response-tap-impulse_concatEp-false', sub); 

get_and_save_ffts(header_avg, data_avg, par.snr_bins, ...
                  'fname', fname, 'save_path', save_path); 

% concatenate epochs
[header_avg, data_avg] = average_trials_per_cond(...
                            header_tap_imp, data_tap_imp, ...
                            'concatenate_epochs', true); 

fname = sprintf('sub-%03d_response-tap-impulse_concatEp-true', sub); 

get_and_save_ffts(header_avg, data_avg, par.snr_bins, ...
                  'fname', fname, 'save_path', save_path); 

                                                                  
