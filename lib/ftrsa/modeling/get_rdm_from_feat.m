function [rdm] = get_rdm_from_feat(feat, varargin)
% Compute RDM from a a set of feature vectors
% 
% Parameters
% ----------
% feat_all_cond : array_like, shape=[n_features, n_conditions] 
%     Feature matrix used to build the RDM. Each column corresponds to the
%     feature vector for one condition. The number of rows corresponds to
%     the number of features. 
% method :  str, {'pearson', 'spearman', 'L2'}, default='pearson'
%     Method that will be used to obtain (dis)similarity of the descriptor
%     vectors across conditions. Either pearson or spearman correlation, or
%     L2 (L2-norm, i.e. Euclidian distance). 
% 
% Returns 
% -------
% rdm : array_like, shape=[n_conditions, n_conditions] 
%     Similarity matrix. 
% 
parser = inputParser; 

addParameter(parser, 'method', 'pearson'); 

parse(parser, varargin{:}); 

method = parser.Results.method; 


% get RDM from the features
if strcmpi(method, 'pearson')
    rdm = corr(feat, 'type', 'Pearson'); 
elseif strcmpi(method, 'spearman')
    rdm = corr(feat, 'type', 'Spearman'); 
elseif strcmpi(method, 'L1')
    distfun_L1 = @(x, y) sum(abs(x - y), 2); 
    rdm = - squareform(pdist(feat', distfun_L1)) + 1; 
elseif strcmpi(method, 'L2')
    rdm = - squareform(pdist(feat')) + 1; 
else
    error('method %s not implementeed', method); 
end

