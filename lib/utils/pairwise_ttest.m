function [p, p_thr, f] = pairwise_ttest(x, varargin)

parser = inputParser(); 

addParameter(parser, 'plot', false); 
addParameter(parser, 'p_adj_method', 'bonf'); 

parse(parser, varargin{:}); 

do_plot = parser.Results.plot; 
p_adj_method = parser.Results.p_adj_method; 

%%

N = size(x, 1); 
n_cond = size(x, 2); 

p = nan(n_cond, n_cond); 

c = 0; 
for i=2:n_cond
    for j=1:i-1
        [~, p(i, j)] = ttest(x(:,i), x(:,j)); 
        c = c+1; 
    end
end

if strcmpi(p_adj_method, 'bonf')
    p = p * c; 
    fprintf('n tests: %d\n', c); 
end

p_thr = p; 
p_thr(p_thr > 0.05) = nan; 


%%

f = []; 
if do_plot
    f = figure; 
    h = imagesc(p < 0.05); 
    set(h, 'AlphaData', ~isnan(p))
    h = pixelgrid; 
    h.Children(1).Color = [0 0 0 ]; 
    h.Children(1).LineWidth = 2;
    h.Children(2).Color = [0 0 0 ]; 
    h.Children(2).LineWidth = 2;
    axis square
    hold on 
    for i=2:n_cond
        for j=1:i-1
            text(j, i, sprintf('%.2g', p(i,j))); 
        end
    end
end
