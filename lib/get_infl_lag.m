function [infl_lags, f1, f2] = get_infl_lag(acf, lags, cycle_dur, varargin)
% temporary function used to find "peaks" in the ACF - atm it's written
% super specifically for 2ioi and 2ioi_stem

parser = inputParser; 

addParameter(parser, 'plot', true); 

parse(parser, varargin{:}); 

do_plot = parser.Results.plot; 

%%

f1 = []; 
f2 = []; 

cols = customcolormap([0, 1], {'#e87b0e', '#6a24a3'}, 13);

infl_lags = nan(1, 13); 

if do_plot
    f1 = figure('Color', 'white', 'Position', [1116 1200 242 1306]); 
    pnl = panel(f1); 
    pnl.pack('v', 13); 
    pnl.de.margin = [0, 1, 0, 0]; 
end

for i=1:13

    idx = dsearchn(lags', cycle_dur * [1/4, 1/2]'); 
    lags_snip = lags([idx(1) : idx(2)]); 
    acf_snip = acf(i, [idx(1) : idx(2)]); 

    d1 = [acf_snip(2:end), nan] - acf_snip; 
    d2 = [d1(2:end), nan] - d1; 

    [~, infl_idx] = findpeaks(-d2, 'npeaks', 1); 

    if isempty(infl_idx)
        infl_idx = dsearchn(lags_snip', cycle_dur/2); 
    end

    infl_lags(i) = lags_snip(infl_idx); 

    if do_plot
        ax = pnl(i).select(); 
        hold on 
        plot([cycle_dur*1/2, cycle_dur*1/2], [-1, 1], 'r--')
        plot([cycle_dur*1/3, cycle_dur*1/3], [-1, 1], 'r--')
        plot([cycle_dur*2/3, cycle_dur*2/3], [-1, 1], 'r--')
        plot(lags, acf(i, :), 'linew', 2, 'color', cols(i,:))
        xlim([0, 0.750])
        ylim([min(acf_snip)-0.1, 1])
        plot(lags_snip(infl_idx), acf_snip(infl_idx), 'ko', 'MarkerFaceColor', 'k'); 
    end

end

infl_lags = (cycle_dur - infl_lags) ./ cycle_dur; 

if do_plot
    f2 = figure('color', 'white', 'Position', [1376 1723 304 251]); 
    plot([1:13], infl_lags, 'k-o', 'linew', 2); 
    xlim([0,14])
    xticks([1, 13])
    ylim([1/2-0.05, 2/3+0.05])
    yticks([1/2, 2/3])
    hold on 
    plot([1,13], [1/2, 1/2], 'r--'); 
    plot([1,13], [2/3, 2/3], 'r--'); 
    box off
end