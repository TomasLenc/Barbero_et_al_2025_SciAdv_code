% This script does all necesary filtering on continuous data, then
% segmentation into trials and downsampling. 

clear 

par = get_par('experiment', '2ioi'); 

sub = 15; 

%%

% load
fpath = fullfile(par.raw_path); 

fname = sprintf('sub-%03d_task-EEG', sub); 

[header, data] = CLW_load(fullfile(fpath, fname)); 

% high-pass filter 
[header, data] = RLW_butterworth_filter(header, data, ...
                                        'filter_type', 'highpass', ...
                                        'low_cutoff', par.filt_cut_low, ...
                                        'filter_order', par.filt_order_low); 

% low-pass filter 
[header, data] = RLW_butterworth_filter(header, data, ...
                                        'filter_type', 'lowpass', ...
                                        'high_cutoff', par.filt_cut_high, ...
                                        'filter_order', par.filt_order_high); 


% segment all trials 
codes = unique({header.events.code}); 

segment_dur = par.trial_dur + 2 * par.epoch_buffer_dur; 

warning('segmenting trials with duration: %.2f s with %.2f s buffer on each side', ...
    par.trial_dur, par.epoch_buffer_dur); 

[header_ep, data_ep] = segment_safe(header, data, codes,...
                                        'x_start', -par.epoch_buffer_dur,...
                                        'x_duration', segment_dur); 

                                    
% save 
fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 
mkdir(fpath); 

header_ep.name =  sprintf('sub-%03d_task-EEG_filt_ep', sub); 
CLW_save(fpath, header_ep, data_ep); 







