% Calcualte ERP waveforms for the EEG response.   

clear

par = get_par('experiment', '2ioi'); 
 
sub = 15; 

% prepare output directories
sub_str = sub_num2str(sub); 

save_path = fullfile(par.deriv_path, 'erp', sub_str); 

mkdir(save_path); 

%% EEG 

% load preprocessed data
fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt*.lw6', sub); 
d = dir(fullfile(fpath, fname)); 
[header, data] = CLW_load(fullfile(d.folder, d.name));

% low-pass filter for ERPs 
[header, data] = RLW_butterworth_filter(header, data, ...
                                      'filter_type', 'lowpass', ...
                                      'high_cutoff', par.filt_erp_cut_low, ...
                                      'filter_order', 4); 

% re-segment 
codes = unique({header.events.code}); 

[header, data] = segment_safe(...
                                header, ...
                                data, ...
                                codes,...
                                'x_start', 0, ...
                                'x_duration', par.trial_dur); 
                            
% rereference 
[header, data] = RLW_rereference(header, data, ...
                                'apply_list', {header.chanlocs.labels}, ...
                                'reference_list', {'TP9', 'TP10'});
      
% get erp chunks
[header_chunk, data_chunk] = chunk_cycles_per_cond(...
                                header, data, par.cycle_dur); 

% decimate the data 
[header_chunk, data_chunk] = RLW_downsample(header_chunk, data_chunk,...
                                      'x_downsample_ratio', par.decim_factor); 
                            
 
% baseline correction (subtract mean of each epoch)
data_chunk = data_chunk - mean(data_chunk, 6); 

% average cycle chunks separtely for each condition 
[header_erp, data_erp] = average_trials_per_cond(...
                            header_chunk, data_chunk); 
                                  
fname = sprintf('sub-%03d_response-erp_lpf-%.0fHz', sub, par.filt_erp_cut_low); 
    
header_erp.name = fname; 

CLW_save(save_path, header_erp, data_erp); 


% % sanity check plot 
% plot_erp_overlay(squeeze(mean(data_erp, 2)), 1/header_erp.xstep)                      

