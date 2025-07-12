function [f, pnl] = plot_mean_bound_boot_distr(bound, ratios_proto_labs, cond_target_proto_idx, max_idx_all_mean_boot)

n_prototypes = length(ratios_proto_labs); 


f = figure('color', 'w', 'pos', [592 524 557 341]); 
pnl = panel(f); 
pnl.pack('h', [10, 70, 10, 10]); 
pnl(2).pack('h', size(max_idx_all_mean_boot, 2)); 
pnl.de.margin = 1; 

hist_edges = linspace(1, n_prototypes, n_prototypes*10); 
for i=1:size(max_idx_all_mean_boot, 2)
    ax = pnl(2, i).select(); 
    h = histogram(ax, max_idx_all_mean_boot(:,i), hist_edges, 'DisplayStyle', 'bar', 'linew', 1);
    h.EdgeColor = 'none'; 
    if i<bound
        h.FaceColor = [135, 175, 201]/255; 
    else
        h.FaceColor = [191, 55, 21]/255; 
    end
    hold on 
    plot(repmat([1:n_prototypes], 2, 1), ...
         [0, max(h.Values)], ...
         'k:', 'linew', 1); 
    plot([cond_target_proto_idx(i), cond_target_proto_idx(i)], ...
         [0, max(h.Values)], ...
         'r:', 'linew', 3); 
     
     
    view(ax, [90 -90]); 
    ax.XLim = [0, n_prototypes+1];
    ax.XTick = []; 
    if i==1
        ax.XTick = [1:n_prototypes]; 
        ax.XTickLabel = ratios_proto_labs; 
    end
    ax.YAxis.Visible = 'off'; 
    ax.XAxis.Direction = 'reverse'; 

end

ax = pnl(1).select(); 
data_for_hist = reshape(max_idx_all_mean_boot(:, 1:bound-1), [], 1); 
h = histogram(ax, data_for_hist, hist_edges, 'DisplayStyle', 'bar', 'linew', 2);
h.EdgeColor = [49, 97, 127] / 255; 
h.FaceColor = [135, 175, 201]/255; 
hold on 
plot(repmat([1:n_prototypes], 2, 1), ...
     [0, max(h.Values)], ...
     'k:', 'linew', 1); 
plot([cond_target_proto_idx(1), cond_target_proto_idx(1)], ...
     [0, max(h.Values)], ...
     'r:', 'linew', 3); 
plot([cond_target_proto_idx(end), cond_target_proto_idx(end)], ...
     [0, max(h.Values)], ...
     'r:', 'linew', 3); 
view(ax, [90 -90]); 
ax.XLim = [0, n_prototypes+1];
ax.XTick = []; 
ax.YAxis.Visible = 'off'; 
ax.XAxis.Direction = 'reverse'; 

ax = pnl(3).select(); 
data_for_hist = reshape(max_idx_all_mean_boot(:, bound:end), [], 1); 
h = histogram(ax, data_for_hist, hist_edges, 'DisplayStyle', 'bar', 'linew', 2);
h.EdgeColor = [168, 11, 0]/255; 
h.FaceColor = [191, 55, 21]/255; 
hold on 
plot(repmat([1:n_prototypes], 2, 1), ...
     [0, max(h.Values)], ...
     'k:', 'linew', 1); 
plot([cond_target_proto_idx(1), cond_target_proto_idx(1)], ...
     [0, max(h.Values)], ...
     'r:', 'linew', 3); 
plot([cond_target_proto_idx(end), cond_target_proto_idx(end)], ...
     [0, max(h.Values)], ...
     'r:', 'linew', 3); 
view(ax, [90 -90]); 
ax.XLim = [0, n_prototypes+1];
ax.XTick = []; 
ax.YAxis.Visible = 'off'; 
ax.XAxis.Direction = 'reverse'; 

ax = pnl(4).select(); 
data_for_hist = reshape(max_idx_all_mean_boot, [], 1); 
h = histogram(ax, data_for_hist, hist_edges, 'DisplayStyle', 'bar', 'linew', 2);
h.EdgeColor = [66, 66, 66]/255; 
h.FaceColor = [204, 204, 204]/255; 
hold on 
plot(repmat([1:n_prototypes], 2, 1), ...
     [0, max(h.Values)], ...
     'k:', 'linew', 1); 
plot([cond_target_proto_idx(1), cond_target_proto_idx(1)], ...
     [0, max(h.Values)], ...
     'r:', 'linew', 3); 
plot([cond_target_proto_idx(end), cond_target_proto_idx(end)], ...
     [0, max(h.Values)], ...
     'r:', 'linew', 3); 
view(ax, [90 -90]); 
ax.XLim = [0, n_prototypes+1];
ax.XTick = []; 
ax.YAxis.Visible = 'off'; 
ax.XAxis.Direction = 'reverse'; 


pnl(1).marginright = 15; 
pnl(2).marginright = 10; 
pnl(3).marginright = 5; 

if n_prototypes > 20
    pnl.fontsize = 6; 
end
