function [x, t] = tap_onset_to_continuous(tap_index, trial_dur, fs, varargin)
% Take a series of tap onset times and generate an onset-only continuous
% signal (remove amplitude related info). Optionally, the signal can be
% convolved with a gaussian kernel for smoothing. 
% 
% Parameters
% ----------
% tap_index : array of integers, shape=[1, n_taps]
%     Indices of the detected tap onsets, each row is one trial. 
% fs : int
%     Sampling rate.
% kernel_fwhm : float, optional, default=[]
%     Full Width at Half Maximum (time-domain) of the smoothing gaussian
%     kernel in seconds. If not provided, no smoothing will be performed.
% 
% Returns 
% -------
% x : shape=[n_trials, time]
%     Input time series. 
%

parser = inputParser; 

addParameter(parser, 'kernel_fwhm', []); 

parse(parser, varargin{:}); 

kernel_fwhm = parser.Results.kernel_fwhm; 

%%

n_trials = size(tap_index, 1); 

N = round(trial_dur * fs); 

x = zeros(n_trials, N); 
for i_trial=1:n_trials
    x(i_trial, tap_index(i_trial, :)) = 1; 
end

% convolve (smooth) with  a gaussian window
if ~isempty(kernel_fwhm)
    x = conv_gauss(x, kernel_fwhm, fs); 
end

% create time vector
t = [0 : N-1] / fs; 