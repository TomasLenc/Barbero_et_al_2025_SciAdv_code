function [tap_onset_times, tap_indices] = extract_taps(x, fs, amplitude_thr, min_asynch, varargin)
% Extract tap times from continuous time series (e.g. recorded using a
% microphone). The extraction is based on thresholding signal's
% amplitude, and minimum inter-tap interval. 
% 
% Parameters
% ----------
% x : shape=[1, time]
%     Input time series. 
% fs : int
%     Sampling rate.
% amplitude_thr : float, optional, default=$
%     Amplitude threshold used to detect taps whenever the signal is above
%     this threshold. 
% min_asynch : float, optional, default=$
%     If the signal amplitude is above `amplitude_thr`, it is only counted
%     as a tap onset if the signal has not been above `amplitude_thr` for
%     at least `min_asynch` seconds before. 
% 
% Returns 
% -------
% tap_onset_times : array of floats 
%     Times in seconds of the detected tap onsets in the input array `x`. 
% tap_index : array of integers 
%     Indices of the detected tap onsets in the input array `x`. 
% 

parser = inputParser; 

addParameter(parser, 'plot_diagnostic', false); 

parse(parser, varargin{:}); 

plot_diagnostic = parser.Results.plot_diagnostic; 

%%

tap_indices = find(x > amplitude_thr) ; 

asy = [Inf, diff(tap_indices) / fs]; 

tap_indices(asy < min_asynch) = []; 

tap_onset_times = tap_indices / fs; 


% sanity check
if plot_diagnostic
    figure('color', 'w', 'pos', [-415 1596 2517 682]); 
    t = [0:length(x)-1]/fs; 
    plot(t, x); 
    hold on 
    plot(t(tap_indices), x(tap_indices), 'ro')   
    plot([0, t(end)], [amplitude_thr, amplitude_thr], 'k:', 'linew', 2)

end

