function [bf10, f] = pairwise_bf_ttest(x, varargin)

parser = inputParser(); 

addParameter(parser, 'plot', false); 

parse(parser, varargin{:}); 

do_plot = parser.Results.plot; 

%%

N = size(x, 1); 
n_cond = size(x, 2); 

bf10 = nan(n_cond, n_cond); 

c = 0; 
for i=2:n_cond
    for j=1:i-1
        bf10(i, j) = bf.ttest(x(:,i), x(:,j));   % paired samples BF
        c = c+1; 
    end
end


%%

f = []; 
if do_plot
    f = figure; 
    h = imagesc((bf10 < 1/3) | (bf10 > 3)); 
    set(h, 'AlphaData', ~isnan(bf10))
    h = pixelgrid; 
    h.Children(1).Color = [0 0 0 ]; 
    h.Children(1).LineWidth = 2;
    h.Children(2).Color = [0 0 0 ]; 
    h.Children(2).LineWidth = 2;
    axis square
    hold on 
    for i=2:n_cond
        for j=1:i-1
            text(j, i, sprintf('%.2g', bf10(i,j))); 
        end
    end
end
