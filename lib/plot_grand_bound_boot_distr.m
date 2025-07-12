function [f, pnl] = plot_grand_bound_boot_distr(bound, ratios_proto_labs, cond_target_proto_idx, max_idx_boot)

% bound = 6; 

n_prototypes = length(ratios_proto_labs); 


f = figure('color', 'w', 'pos', [951 604 382 303]); 

if n_prototypes > 20
   f.Position(4) = 500;  
end

pnl = panel(f); 
pnl.pack('h', [10, 70, 10, 10]); 
pnl(2).pack('h', size(max_idx_boot, 2)); 
pnl.de.margin = 1; 

hist_edges = [0.5 : 1 : n_prototypes+0.5]; 
for i=1:13
    ax = pnl(2, i).select(); 

    
    h = histogram(ax, max_idx_boot(:,i), hist_edges, 'DisplayStyle', 'bar', 'linew', 2);
     
    hold(ax, 'on'); 
    plot(repmat([1:n_prototypes], 2, 1), ...
         [0, max(h.Values)], ...
         'k:', 'linew', 0.1); 
    plot([cond_target_proto_idx(i), cond_target_proto_idx(i)], ...
         [0, max(h.Values)], ...
         '-', 'linew', 2, 'color', [112, 112, 112]/255); 
    delete(h'); 
    h = histogram(ax, max_idx_boot(:,i), hist_edges, 'DisplayStyle', 'bar', 'linew', 2);
    h.EdgeColor = 'none'; 
    if i<=bound
        h.FaceColor = [135, 175, 201]/255; 
    else
        h.FaceColor = [191, 55, 21]/255; 
    end
     
    view(ax, [90 -90]); 
    ax.XLim = [0, n_prototypes+1];
    ax.XTick = []; 
    if i==1
        ax.XTick = [1:n_prototypes]; 
        ax.XTickLabel = ratios_proto_labs; 
        ax.TickLength = [0,0];  
    end
    ax.YAxis.Visible = 'off'; 
    ax.XAxis.Direction = 'reverse'; 
end

ax = pnl(1).select(); 
data_for_hist = reshape(max_idx_boot(:, 1:bound), [], 1); 
h = histogram(ax, data_for_hist, hist_edges, 'DisplayStyle', 'bar', 'linew', 2);
h.EdgeColor = [49, 97, 127] / 255; 
h.FaceColor = [135, 175, 201]/255; 
view(ax, [90 -90]); 
ax.XLim = [0, n_prototypes+1];
ax.XTick = []; 
ax.YAxis.Visible = 'off'; 
ax.XAxis.Direction = 'reverse'; 

ax = pnl(3).select(); 
data_for_hist = reshape(max_idx_boot(:, bound+1:end), [], 1); 
h = histogram(ax, data_for_hist, hist_edges, 'DisplayStyle', 'bar', 'linew', 2);
h.EdgeColor = [168, 11, 0]/255; 
h.FaceColor = [191, 55, 21]/255; 
view(ax, [90 -90]); 
ax.XLim = [0, n_prototypes+1];
ax.XTick = []; 
ax.YAxis.Visible = 'off'; 
ax.XAxis.Direction = 'reverse'; 

% ax = pnl(4).select(); 
% data_for_hist = reshape(max_idx_boot, [], 1); 
% h = histogram(ax, data_for_hist, hist_edges, 'DisplayStyle', 'bar', 'linew', 2);
% h.EdgeColor = [66, 66, 66]/255; 
% h.FaceColor = [204, 204, 204]/255; 
% view(ax, [90 -90]); 
% ax.XLim = [0, n_prototypes+1];
% ax.XTick = []; 
% ax.YAxis.Visible = 'off'; 
% ax.XAxis.Direction = 'reverse'; 


pnl(1).marginright = 10; 
pnl(2).marginright = 4; 
pnl(3).marginright = 4; 

if n_prototypes > 20
    pnl.fontsize = 6; 
end