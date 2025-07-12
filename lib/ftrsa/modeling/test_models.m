function [r, p, p_adj] = test_models(rdm, models, varargin)
% Compute correlation of a response RDM with a model RDM. In principle, any
% RDM can be passed as a model (e.g. if we want to correlate neural and
% behaviroal RDMs). Spearman correlation is always used. 
% 
% Parameters
% ----------
% rdm : array_like, shape=[n, n]
%     Representational dissimiliary matrix with n dimensions. 
% models : (cell array of) array_like, shape=[n, n]
%     Cell array containing all models that are going to be tested. The
%     number of models must be equal to the number of rows in `models_info`
%     table. Can be also a single model passed as an [n x n] numeric array.
% covariate :  array_like, shape=[n, n], optional
%     DSM that will be regressed out before correlating rdm and categorical 
%     model. 
% 
% Returns 
% -------
% r : array_like, shape=[n_models, 1]
%     List of correlation coefficients, separately for each model. 
% p : array_like, shape=[n_models, 1]
%     List of p-values corresponding to the r values, separately for each
%     model.
% p_adj : array_like, shape=[n_models, 1]
%     List of p-values corresponding to the r values, separately for each
%     model, adjusted for the number of tested models (Bonferroni
%     correction, i.e., the original pvalue is multiplied by the number of
%     tested models). 
% 

parser = inputParser; 

addParameter(parser, 'covariate', []); 
addParameter(parser, 'verbose', true); 

parse(parser, varargin{:}); 

% check if we got multiple models in a cell array or a single model as an
% array
if ~iscell(models)
    assert(isequal(ndims(models), 2), ...
        'Incorrect model dimensions. The model must be a NxN matrix!')
    assert(size(models, 1) == size(models, 2), ...
        'Incorrect model dimensions. The model must be a NxN matrix!')
    models = {models}; 
    
end

n_models = length(models); 

r = nan(n_models, 1); 
p = nan(n_models, 1); 

for i_model=1:length(models)
    
    model = models{i_model}; 

    [r(i_model), p(i_model)] = get_rdm_model_corr(rdm, model, varargin{:}); 
    
end

% adjust pvalue
p_adj = min(p * n_models, 1); 



