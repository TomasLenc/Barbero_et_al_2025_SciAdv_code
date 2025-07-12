function [r, p] = get_rdm_model_corr(rdm, model, varargin)
% Compute correlation of a response RDM with a model RDM. In principle, any
% RDM can be passed as a model (e.g. if we want to correlate neural and
% behaviroal RDMs). Spearman correlation is always used. 
% 
% Parameters
% ----------
% rdm : array_like, shape=[n, n]
%     Representational dissimiliary matrix with n dimensions. 
% model : array_like, shape=[n, n]
%     RDM of the model to evaluate. 
% covariate :  array_like, shape=[n, n], optional
%     DSM that will be regressed out before correlating rdm and categorical 
%     model. 
% 
% Returns 
% -------
% r : float 
%     Correlation between the rdm and model. 
% p : float 
%     Pvalue for the correlation. 
% 

parser = inputParser; 

addParameter(parser, 'covariate', []); 
addParameter(parser, 'verbose', 0); 

parse(parser, varargin{:}); 

covariate_rdm = parser.Results.covariate; 
verbose = parser.Results.verbose; 

%----------------------------------------------------------------------
if verbose && ~isempty(covariate_rdm)
    warning('running partial correlation (regressing out covariate RDM)'); 
end

n = size(rdm, 1); 
mask = tril(ones(n, n, 'logical')) & ~diag(ones(1, n, 'logical')) ; 

rdm_lower_tri = rdm(mask); 
model_lower_tri = model(mask); 

if ~isempty(covariate_rdm)
    % run partial correlation
    covariate_lower_tri = covariate_rdm(mask); 

    [r, p] = partialcorr(...
                        model_lower_tri,...
                        rdm_lower_tri, ...
                        covariate_lower_tri, ...
                        'Type','Spearman'); 
else
    % run correlation
    [r, p] = corr(...
                model_lower_tri,...
                rdm_lower_tri, ...
                'Type','Spearman'); 
end


