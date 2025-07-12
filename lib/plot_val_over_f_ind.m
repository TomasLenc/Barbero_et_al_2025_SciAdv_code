function [f, ax] = plot_val_over_f_ind(freq, val, varargin)

parser = inputParser; 

addParameter(parser, 'ax', []); 
addParameter(parser, 'figsize', [809 666 105 203]); 
addParameter(parser, 'p', []); 
addParameter(parser, 'add_line', false); 
addParameter(parser, 'color', [0, 0, 0] / 255); 

parse(parser, varargin{:}); 

ax = parser.Results.ax; 
c = parser.Results.color; 
add_line = parser.Results.add_line; 
figsize = parser.Results.figsize; 
p = parser.Results.p; 

f = []; 
if isempty(ax)
    f = figure('pos', figsize, 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select(); 
end

r_mu = mean(val, 1); 
r_sem = std(val, [], 1) / sqrt(size(val, 1)); 

linespec = 'o'; 
if add_line
    linespec = '-o'; 
end

plot(freq, r_mu, linespec, 'linew', 2, 'color', c, 'markerfacecolor', c, 'markersize', 5)
hold on 

for i=1:length(r_mu)
    
    plot([freq(i), freq(i)], ...
         [r_mu(i)-r_sem(i), r_mu(i)+r_sem(i)], ...
          'color', c, 'linew', 2)
    
end

if ~isempty(p)
    mask_signif = p < 0.05; 
    plot(freq(mask_signif), r_mu(mask_signif) + r_sem(mask_signif) + 0.02, 'k*');  
end
