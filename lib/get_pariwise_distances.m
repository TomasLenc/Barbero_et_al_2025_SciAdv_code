function [dist_matrix, f, ax, max_idx] = get_pariwise_distances(proto, resp, varargin)

parser = inputParser; 

addParameter(parser, 'plot', false)
addParameter(parser, 'plot_max', false)
addParameter(parser, 'normalize_columns', false)
addParameter(parser, 'ratios_proto_labs', [])

parse(parser, varargin{:}); 

do_plot = parser.Results.plot; 
plot_max = parser.Results.plot_max; 
normalize_columns = parser.Results.normalize_columns; 
ratios_proto_labs = parser.Results.ratios_proto_labs; 


n_cond_proto = size(proto, 2); 
n_cond_resp = size(resp, 2); 

dist_matrix = nan(n_cond_proto, n_cond_resp); 

for i_cond_proto=1:n_cond_proto
    for i_cond_resp=1:n_cond_resp
        dist_matrix(i_cond_proto, i_cond_resp) = ...
                    corr(proto(:, i_cond_proto), ...
                         resp(:, i_cond_resp), ...
                         'type', 'pearson'); 
    end
end

if normalize_columns
    % re-normalize each column of the matrix to max 1
    dist_matrix = dist_matrix ./ max(abs(dist_matrix), [], 1); 
end

% find maximum in each column 
[~, max_idx] = max(dist_matrix, [], 1); 

if do_plot
    fig_size = [680 722 229 256]; 
    f = figure('color', 'w', 'Position', fig_size); 
    ax = axes(f); 
    plot_rdm(dist_matrix, 'ax', ax); 
    xlabel(ax, 'response'); 
    ylabel(ax, 'prototype'); 
    ax.XAxis.Visible = 'on'; 
    ax.YAxis.Visible = 'on'; 
    ax.CLim = [-max(abs(dist_matrix(:))), max(abs(dist_matrix(:)))]; 
    colormap(brewermap(256, '-RdYlBu'))
    
    if plot_max
        if n_cond_proto > 20
            marker_size = 4; 
        else
            marker_size = 8; 
        end
        plot([1:size(dist_matrix,2)], max_idx, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', marker_size)
    end

    if ~isempty(ratios_proto_labs)
        ax.YTick = [1 : length(ratios_proto_labs)]; 
        ax.YTickLabel = ratios_proto_labs; 
    end
else
    f = []; 
    ax = [];    
end

if n_cond_proto > 20
    ax.FontSize = 6; 
end

if n_cond_proto > 40
    ax.FontSize = 3; 
end