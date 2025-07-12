function [f, pnl] = plot_mX_lw(header_mX, data_mX, chan_idx, frex, varargin)
% Plot magnitude spectrum for each condition, using lw6 data input. 
% 
% Parameters
% ----------
% header_mX : struct
%     Lw6 header. 
% data_mX : array_like, shape=[conditions, channels, 1, 1, 1, frequency]
%     Lw6 frequency domain data with magnitude spectra. 
% chan_idx : array_like
%     Integer indices of the channels that will be used in the plot.
%     Magnitude spectra across these channels will be averaged for each
%     condition and plotted. 
% frex : array_like
%     Frequencies of interests. These will be highlighted in the plot. 
% z_snr :  array_like, optional
%     SNR value for each condition. If passed, it will be plotted as text
%     label above each spectrum. 
% fmin :  float, optional, default=0.3
%     Minimum frequency limit of the plot. 
% fmax :  float, optional, default=max(frex)+1
%     Maximum frequency to plot. 
% 
% Returns 
% -------
% f : figure object 
%     Generated figure. 

parser = inputParser(); 

addParameter(parser, 'pnl', []); 
addParameter(parser, 'z_snr', []); 
addParameter(parser, 'fmin', 0.3); 
addParameter(parser, 'fmax', max(frex) + 1); 

parse(parser, varargin{:}); 

pnl = parser.Results.pnl; 
z_snr = parser.Results.z_snr; 
fmin = parser.Results.fmin; 
fmax = parser.Results.fmax; 


n_cond = size(data_mX, 1);

mX = cell(n_cond, 1); 

for i_cond=1:n_cond
                
    frex_idx = round((frex - header_mX.xstart) / header_mX.xstep) + 1; 
    min_freq_idx = round((fmin - header_mX.xstart) / header_mX.xstep) + 1; 
    max_freq_idx = round((fmax - header_mX.xstart) / header_mX.xstep) + 1; 
        
    mX{i_cond} = ensure_row(squeeze(mean(data_mX(i_cond, chan_idx, 1, 1, 1, :), 2))); 
    
end

freq = [0 : header_mX.datasize(end)-1] * header_mX.xstep + header_mX.xstart;

[f, pnl] = plot_mX(mX, freq, frex, varargin{:}); 

tit = sprintf('channels: %s', ...
    chan_list_to_str({header_mX.chanlocs(chan_idx).labels})); 

pnl.title(tit); 


