function [f, pnl] = plot_grand_bound_boot_distr_h(bound, ratios_proto, ratios_stim, max_idx_boot, varargin)

parser = inputParser; 
addParameter(parser, 'show_proto_ticks', false);
addParameter(parser, 'show_stim_ticks', false);
addParameter(parser, 'highlight_ratios', []);

parse(parser, varargin{:});

show_proto_ticks = parser.Results.show_proto_ticks;
show_stim_ticks = parser.Results.show_stim_ticks;
highlight_ratios = parser.Results.highlight_ratios;


%%

x_kde = [0 : 0.001 : 1]; 

lighten_color = @(col, factor) col + (1 - col) * factor;

%%

f = figure('color', 'w', 'pos', [571 365 270 452]); 

pnl = panel(f); 
pnl.pack('v', [15, 70, 15]); 
pnl(2).pack('v', size(max_idx_boot, 2)); 
pnl.de.margin = 1; 

% convert from index to ratio 
max_ratio_boot = ratios_proto(max_idx_boot); 

proto_ratio_step = ratios_proto(2) - ratios_proto(1); 

hist_edges = [...
    ratios_proto(1)-proto_ratio_step/2 : ...
    proto_ratio_step : ...
    ratios_proto(end)+proto_ratio_step/2 ...
    ]; 

n_prototypes = length(ratios_proto); 

cond_target_proto_idx = dsearchn(ensure_col(ratios_proto), ensure_col(ratios_stim)); 

for i_cond=1:13
    
    ax = pnl(2, i_cond).select(); 

    h = histogram(ax, max_ratio_boot(:,i_cond), hist_edges, ...
                  'Normalization', 'pdf', ...
                  'DisplayStyle', 'bar', ...
                  'linew', 2);
     
    hold(ax, 'on'); 
    
    if show_proto_ticks       
        plot(repmat(ensure_row(ratios_proto), 2, 1), ...
             [0, max(h.Values)], ...
             'k:', 'linew', 0.1); 
    end
    if show_stim_ticks       
         plot([ratios_stim(i_cond), ratios_stim(i_cond)], ...
             [0, max(h.Values)], ...
             '-', 'linew', 5, 'color', [112, 112, 112]/255); 
    end
    
    delete(h'); 
    h = histogram(ax, max_ratio_boot(:,i_cond), hist_edges, ...
                  'Normalization', 'pdf', ...
                  'DisplayStyle', 'bar', ...
                  'linew', 2);
              
    h.EdgeColor = 'none'; 
    
    if i_cond<=bound
        col_light = [135, 175, 201]/255; 
        col_dark = [49, 97, 127] / 255; 
    else
        col_light = [191, 55, 21]/255; 
        col_dark = [168, 11, 0]/255; 
    end
    
    h.FaceColor = lighten_color(col_light, 0.5); 

    ax.YAxis.Visible = 'off'; 
    ax.XLim = [hist_edges(1), hist_edges(end)];
    ax.XTick = []; 
    
    % KDE     
    kde = ksdensity(max_ratio_boot(:,i_cond), x_kde, 'bandwidth', 0.005); 
    plot(x_kde, kde, 'color', lighten_color(col_dark, 0.5), 'linew', 2)
        
end


%% marginal for caterogy 1

col_light = [135, 175, 201]/255; 
col_dark = [49, 97, 127] / 255; 

ax = pnl(1).select(); 
data_for_hist = reshape(max_ratio_boot(:, 1:bound), [], 1); 
h = histogram(ax, data_for_hist, hist_edges, ...
              'Normalization', 'pdf', ...
              'DisplayStyle', 'bar', ...
              'linew', 2);
h.EdgeColor = 'none'; 
h.FaceColor = col_light; 
ax.XLim = [hist_edges(1), hist_edges(end)];ax.XTick = []; 
ax.YAxis.Visible = 'off'; 

hold(ax, 'on'); 
if ~isempty(highlight_ratios)
    for i_cond=1:length(highlight_ratios)
       plot(ax, [highlight_ratios(i_cond), highlight_ratios(i_cond)], ...
                [0, max(h.Values)*1.1], 'm', 'linew', 4); 
    end
end

% KDE     
kde = ksdensity(data_for_hist, x_kde, 'bandwidth', 0.005); 
plot(x_kde, kde, 'color', col_dark, 'linew', 2)

%% marginal category 2

col_light = [191, 55, 21]/255; 
col_dark = [168, 11, 0]/255; 

ax = pnl(3).select(); 
data_for_hist = reshape(max_ratio_boot(:, bound+1:end), [], 1); 
h = histogram(ax, data_for_hist, hist_edges, ...
              'Normalization', 'pdf', ...
              'DisplayStyle', 'bar', ...
              'linew', 2);
h.EdgeColor = 'none'; 
h.FaceColor = col_light; 
ax.XLim = [hist_edges(1), hist_edges(end)];
ax.XTick = []; 
ax.YAxis.Visible = 'off'; 

hold(ax, 'on'); 
if ~isempty(highlight_ratios)
    for i_cond=1:length(highlight_ratios)
       plot(ax, [highlight_ratios(i_cond), highlight_ratios(i_cond)], ...
                [0, max(h.Values)*1.1], 'm', 'linew', 4); 
    end
end

% KDE     
kde = ksdensity(data_for_hist, x_kde, 'bandwidth', 0.005); 
plot(x_kde, kde, 'color', col_dark, 'linew', 2)


%%

ax.XTick = ratios_proto; 
ax.XTickLabel = get_ratio_labs(ratios_proto, 'n_decim', 3, 'no_zero', true); 
ax.TickLength = [0,0];  
ax.XTickLabelRotation = -60; 

pnl(1).marginbottom = 4; 
pnl(2).marginbottom = 4; 
pnl(3).marginbottom = 4; 

pnl.margin = [5, 5, 5, 5]; 

if n_prototypes > 20
    pnl.fontsize = 4; 
end