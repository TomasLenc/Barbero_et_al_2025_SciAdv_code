function [rdm, f] = rsa_stim_erp(varargin)

parser = inputParser(); 

addParameter(parser, 'plot', false); 
addParameter(parser, 'save', false); 
addParameter(parser, 'par', []); 
addParameter(parser, 'lfp_cutoff', []); 
addParameter(parser, 'flip_rdm', true); 

parse(parser, varargin{:}); 

do_plot = parser.Results.plot; 
do_save = parser.Results.save; 
par = parser.Results.par; 
lfp_cutoff = parser.Results.lfp_cutoff; 
flip_rdm = parser.Results.flip_rdm; 

if isempty(par)
    par = get_par(); 
end

if isempty(lfp_cutoff)
    lfp_cutoff = par.filt_erp_cut_low; 
end

% load
fpath = fullfile(par.deriv_path, 'stimuli'); 
load(fullfile(fpath, 'stim.mat')); 

% envelope
stim_env = cellfun(@(x) abs(hilbert(x)), stim.s, 'uni', 0); 

% low-pass filter 
[b, a] = butter(4, lfp_cutoff / (stim.fs / 2), 'low'); 
stim_env = cellfun(@(x) filtfilt(b,a,x), stim_env, 'uni', 0); 

% decimate
decim_factor = 100; 
fs = stim.fs / decim_factor; 
stim_env = cellfun(@(x) x(1:decim_factor:end), stim_env, 'uni', 0); 

% chunk
stim_env_chunk = cellfun(@(x) epoch_chunks(x, fs, par.cycle_dur), stim_env, 'uni', 0); 
stim_env_erp = cellfun(@(x) mean(x, 1), stim_env_chunk, 'uni', 0); 

stim_env_erp = cat(1, stim_env_erp{:}); 

rdm = get_rdm_from_feat(stim_env_erp', 'method', 'pearson'); 

if do_save 
    save_path = fullfile(par.deriv_path, 'results', ...
                         'response-stim-erp'); 

    mkdir(save_path); 
    
    fname = sprintf('response-stim-env_feat-erp_method-pearson_rdm_lpf-%.0fHz',...
                    lfp_cutoff); 
    save(fullfile(save_path, [fname, '.mat']), 'rdm'); 
end
  
f = []; 
if do_plot
    f = plot_rdm(rdm, 'ratios', stim.ratios, 'flip', flip_rdm); 
    if do_save 
        print(f, '-dpng', '-painters', '-r300', fullfile(save_path, fname))
    end
end






