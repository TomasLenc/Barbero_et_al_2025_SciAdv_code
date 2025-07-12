function [f, pnl] = plot_dist_matrix_distr(bound, ratios_proto_labs, cond_target_proto_idx, dist_matrix)

n_prototypes = length(ratios_proto_labs); 

f = figure('color', 'w', 'pos', [951 604 382 303]); 

if n_prototypes > 20
   f.Position(4) = 500;  
end

% rectify to remove negative correlations
dist_matrix = max(dist_matrix, 0); 

pnl = panel(f); 
pnl.pack('h', [10, 70, 10, 10]); 
pnl(2).pack('h', size(dist_matrix, 2)); 
pnl.de.margin = 1; 

for i=1:13
    ax = pnl(2, i).select(); 
  
    hold(ax, 'on'); 
    plot(repmat([1:n_prototypes], 2, 1), ...
         [min(dist_matrix(:,i)), max(dist_matrix(:,i))], ...
         'k:', 'linew', 0.1); 
    plot([cond_target_proto_idx(i), cond_target_proto_idx(i)], ...
         [min(dist_matrix(:,i)), max(dist_matrix(:,i))], ...
         '-', 'linew', 2, 'color', [112, 112, 112]/255); 
     
    h = bar(ax, dist_matrix(:,i), 'LineStyle', 'none');
     
    if i<=bound
        h.FaceColor = [135, 175, 201]/255; 
    else
        h.FaceColor = [191, 55, 21]/255; 
    end
    h.FaceAlpha = 0.5; 
     
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
data = mean(dist_matrix(:, 1:bound), 2); 
h = bar(ax, data, 'LineStyle', 'none');
h.EdgeColor = [49, 97, 127] / 255; 
h.FaceColor = [135, 175, 201]/255; 
view(ax, [90 -90]); 
ax.XLim = [0, n_prototypes+1];
ax.XTick = []; 
ax.YAxis.Visible = 'off'; 
ax.XAxis.Direction = 'reverse'; 

ax = pnl(3).select(); 
data = mean(dist_matrix(:, bound+1:end), 2); 
h = bar(ax, data, 'LineStyle', 'none');
h.EdgeColor = [168, 11, 0]/255; 
h.FaceColor = [191, 55, 21]/255; 
view(ax, [90 -90]); 
ax.XLim = [0, n_prototypes+1];
ax.XTick = []; 
ax.YAxis.Visible = 'off'; 
ax.XAxis.Direction = 'reverse'; 


pnl(1).marginright = 10; 
pnl(2).marginright = 4; 
pnl(3).marginright = 4; 

if n_prototypes > 20
    pnl.fontsize = 6; 
end