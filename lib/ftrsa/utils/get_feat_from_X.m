function [feat_all_cond] = get_feat_from_X(X, freq, frex, varargin)
% Compute feature vectors from complex DFT. The input is a complex DFT, and
% the function extracts real and imaginary values at a set of frequenies
% (frex), concatenates them into long feature vectors (also adding any
% other dimensions that are in the input DFT array, e.g. channels. The
% output is a matrix where each column is a feature vector for one
% condition.
% 
% Parameters
% ----------
% X : cell array | array_like
%     All the data. If cell array is passed, each condition corresponds to
%     one cell in the array, and contains the complex fourier spectrum with
%     frequency as the last dimensions. Any dimensions that are before
%     correspond to additional features that will be taken to built the
%     feature vectors and the RDM. If a numeric array is passed, then
%     conditions must be along the first dimension, and frequency along the
%     last dimension. Again, any dimensions in between will be treated as
%     additional features (e.g. EEG channels), and taken separately to
%     build the feature vector for each condition. 
% freq :  array_like
%     Frequency vector. Must have the same length as the last dimension of
%     the DFT array for each condition. 
% frex :  array_like
%     Frequencies of interest in Hz. The real and imaginary components at
%     these frequencies will be used to build the descriptor vector for
%     each condition.
% f_center : float
%     If non-zero, the frequencies in `frex` will be extracted from
%     symmetrical sidebands around `f_center`. For example, if `f_center`
%     is 100, and `frex` is [1, 2, 3] Hz, then the RDM will be built from
%     descriptor vectors corresponding to real and imaginary components
%     taken from [97, 98, 99, 101, 102, 103] Hz. 
% feat : str, {'X', 'mX'}, optional, default='X'
%     Name of the feature that will be extradted from the complex spectrum.
%     Default is 'X', which will make feature vectors by concatenating the
%     real and imaginary component of the complex coefficient at each frex.
%     If 'mX', the magnitude (absolute value) will be used instead of the
%     real and imag representation. 
% 
% Returns 
% -------
% feat_all_cond : array_like, shape=[n_conditions, n_conditions] 
%     Feature matrix used to build the RDM. Each column corresponds to the
%     feature vector for one condition. The number of rows corresponds to
%     the number of features. 
%

parser = inputParser; 

addParameter(parser, 'f_center', 0); 
addParameter(parser, 'feat', 'X'); 

parse(parser, varargin{:}); 

f_center = parser.Results.f_center; 
feat_name = parser.Results.feat; 


%%

if f_center ~= 0
    frex = f_center + [-flip(frex), frex]; 
end
    
frex_idx = dsearchn(ensure_col(freq), ensure_col(frex))'; 

% If X is a cell, then each element should be an array containing complex
% DFT spectrum, with frequency at the last dimension. The shape of each
% array should be exactly the same! 
if iscell(X)
    % each cell element is an array that can have many dimensions (e.g.
    % channels), but the last one must be frequency. 
    assert(length(unique(cellfun(@(x) numel(x), X, 'uni', 1))) == 1, ...
        'every condition in X must be an array with the same size!'); 
    
    X = cat(1, X{:});           
end

% sanity check that the input is complex spectrum
if isreal(X) && ~strcmp(feat_name, 'mX')
    warn_state = warning; 
    warning('on'); 
    warning('Seems like you provided magnitude spectrum, but we need complex! Make sure you know what you are doing!'); 
    warning(warn_state.state); 
end
    
% get features from the DFT. If there are additional dimensions 
n_cond = size(X, 1); 
shape = size(X); 

% Check if there are any dimensions between the first (conditions) and last
% (frequency). If so, we will need bigger feature vectors. 
if ndims(X) > 2 && any(shape(2:end-1) > 1)
    other_dims_shape = shape(2:end-1); 
    n_other_dims = sum(other_dims_shape(other_dims_shape > 1)); 
else
    n_other_dims = 1; 
end

if strcmp(feat_name, 'X')
    feat_all_cond = nan(n_other_dims * 2 * length(frex), n_cond); 
elseif strcmp(feat_name, 'mX')
    feat_all_cond = nan(n_other_dims * length(frex), n_cond); 
else
    error('feature name not recognized'); 
end


% Go over conditions and extract features at frex. 
for i_cond=1:n_cond
    index = repmat({':'}, 1, ndims(X)); 
    index{1} = i_cond; 
    index{end} = frex_idx; 
    feat = X(index{:}); 
    feat = reshape(feat, [], 1); 
    if strcmp(feat_name, 'X')
        feat_all_cond(:, i_cond) = [real(feat); imag(feat)]; 
    elseif strcmp(feat_name, 'mX')
        feat_all_cond(:, i_cond) = abs(feat); 
    end
end

