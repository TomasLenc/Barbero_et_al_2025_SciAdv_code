% Calcualte "ERP" waveforms for the tapping response.   

clear

par = get_par('experiment', '2ioi'); 

sub = 15; 

%%

% prepare output directories
sub_str = sub_num2str(sub); 

save_path = fullfile(par.deriv_path, 'erp', sub_str); 
mkdir(save_path); 

% load data 
if strcmp(par.experiment, '2ioi')
    tap_response_type = 'force'; 
    decim_factor = 4; 
else
    tap_response_type = 'env'; 
    decim_factor = 100; 
end

fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
fname =  sprintf('sub-%03d_task-TAP_response-%s.lw6', sub, tap_response_type); 
d = dir(fullfile(fpath, fname)); 

[header, data] = CLW_load(fullfile(d.folder, d.name));

% low-pass filter     
[header, data] = RLW_butterworth_filter(header, data, ...
                                                  'filter_type', 'lowpass', ...
                                                  'high_cutoff', par.filt_tap_erp_cut_low, ...
                                                  'filter_order', 4); 

% segment all trials
codes = unique({header.events.code}); 
[header, data] = segment_safe(...
                                header, ...
                                data, ...
                                codes,...
                                'x_start', 0, 'x_duration', ...
                                par.trial_dur); 

% segment into cycle-long chunks                                              
[header_chunk, data_chunk] = chunk_cycles_per_cond(...
                                header, data, par.cycle_dur); 

% decimate the data 
[header_chunk, data_chunk] = RLW_downsample(header_chunk, data_chunk,...
                                      'x_downsample_ratio', decim_factor); 
                            
 
% baseline correction (subtract mean of each epoch)
data_chunk = data_chunk - mean(data_chunk, 6); 

% average cycle chunks separtely for each condition 
[header_erp, data_erp] = average_trials_per_cond(...
                            header_chunk, data_chunk); 
                                  
fname = sprintf('sub-%03d_response-tap-force-erp_lpf-%.0fHz', sub, par.filt_tap_erp_cut_low); 
    
header_erp.name = fname; 

CLW_save(save_path, header_erp, data_erp); 

% % sanity check 
% plot_erp_overlay(squeeze(mean(data_erp, 2)), 1/header_erp.xstep)                      

