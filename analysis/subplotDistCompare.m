function subplotDiscCompare(nrows, ncols, row, plotInfo, gstats, fieldname, imidx)

    barEdgeColor = 'none';
    margin = 0.08;

    nPlots = numel(fieldname);

    if nargin == 6
        imidx = 1;
    end
    
    nBins = length(gstats{imidx}.(fieldname{1}).bins);
    dist = zeros([nPlots nBins]);
    bins = zeros([nPlots nBins]);
    
    for i = 1:nPlots
        dist(i,:) = gstats{imidx}.(fieldname{i}).dist;
        bins(i,:) = gstats{imidx}.(fieldname{i}).bins;
    end
    
    % distribution
    %-------------------------
    
    subplot_tight(nrows, ncols, 2*(row-1) + 1, margin)

    bar(bins', dist');
    xi = min(bins(1,:)); %xi = min(X{1});
    xf = max(bins(1,:)); %xf = max(X{1});
    axis([xi xf 0 1.2*max(dist(:))]);
    
    legend(fieldname, 'Interpreter', 'none');
    if isfield(plotInfo, 'xlabel') xlabel(plotInfo.xlabel); end
    ylabel([]);
    
    title([plotInfo.plotTitle ' distribution'], 'Interpreter', 'tex', 'FontSize', 12);
    

    % cumulative distribution
    %-------------------------
    
    cmap = lines(nPlots);
    subplot_tight(nrows, ncols, 2*row, margin)

    plot(bins(i,:), cumsum(dist(1,:)), 'Color', cmap(1,:), 'LineWidth', 2)
    hold on;
    for i = 2:nPlots
        plot(bins(i,:), cumsum(dist(i,:)), 'Color', cmap(i,:), 'LineWidth', 2);
    end
    hold off;
    axis([xi xf 0 1]);
    if isfield(plotInfo, 'xlabel') xlabel(plotInfo.xlabel); end
    legend(fieldname, 'Location', 'SouthEast', 'Interpreter', 'none');
    title([plotInfo.plotTitle ' cumulative distr'], 'Interpreter', 'tex', 'FontSize', 12);

end