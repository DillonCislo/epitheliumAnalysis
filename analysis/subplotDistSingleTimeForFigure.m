function subplotDistSingleTimeForFigure(nrows, ncols, row, plotInfo, gstats, fieldname, names)
    % call with:
    % subplotDistSingleTimeForFigure(nrows, ncols, row, plotInfo,...
    %                         gstats(INDEX), fieldname, names(INDEX))

    nImages = numel(names);
    
	% are datasets batched? if not make single batch
	if nargin == 7
        batchIdx = {1:nImages};
        batchLabel = {};
        nBatches=1;
        goodseg = 1:nImages;
    else
        nBatches = numel(batchLabel);
        goodseg = cat(2,batchIdx{:});
    end
                    
    barEdgeColor = 'none';
    margin = 0.08;

    nBins = length(gstats{1}.(fieldname).bins);
    dist = zeros([nImages nBins]);
    bins = zeros([nImages nBins]);
    
    for i = goodseg
        dist(i,:) = gstats{i}.(fieldname).dist;
        bins(i,:) = gstats{i}.(fieldname).bins;
    end
    
    if isempty(batchLabel)
        cmap = jet(nImages);
    else
        subcmap = jet(nBatches);
        cmap = zeros([nImages 3]);
        for i = 1:nBatches
            for j = batchIdx{i}
                cmap(j,:) = subcmap(i,:);
            end
        end
    end
    
    % distribution
    %-------------------------
    
    subplot_tight(nrows, ncols, 2*(row-1) + 1, margin)

    hold on;
    n = 1;
    for i = 1:nBatches
        for j = 1:numel(batchIdx{i})
            h(n) = plot(bins(n,:), dist(n,:), 'Color', 'k', 'LineWidth', 2);
            n = n+1;
        end
    end
    hold off;
    
    xi = min(bins(1,:)); %xi = min(X{1});
    xf = max(bins(1,:)); %xf = max(X{1});
    axis([xi xf 0 1.2*max(dist(:))]);
    
%     if isempty(batchLabel)
%         legend(names, 'Interpreter', 'none');
%     else
%         reps = [];
%         for i = 1:nBatches
%             if ~isempty(batchIdx{i})
%                 reps = [reps h(batchIdx{i}(1))];
%             end
%         end
%         legend(reps, batchLabel);
%     end
    if isfield(plotInfo, 'xlabel') xlabel(plotInfo.xlabel); end
    ylabel([]);
    
    title([plotInfo.plotTitle ' distribution'], 'Interpreter', 'tex', 'FontSize', 12);
    

    % cumulative distribution
    %-------------------------
    
    subplot_tight(nrows, ncols, 2*row, margin)
    
    hold on;
    h = [];
    n = 1;
    for i = 1:nBatches
        for j = 1:numel(batchIdx{i})
            h(n) = plot(bins(n,:), cumsum(dist(n,:)), 'Color', 'k', 'LineWidth', 2);
            n = n+1;
        end
    end
    hold off;
    
    axis([xi xf 0 1]);
    if isfield(plotInfo, 'xlabel') xlabel(plotInfo.xlabel); end
%     if isempty(batchLabel)
%         legend(names, 'Location', 'SouthEast', 'Interpreter', 'none');
%     else
%         reps = [];
%         for i = 1:nBatches
%             if ~isempty(batchIdx{i})
%                 reps = [reps h(batchIdx{i}(1))];
%             end
%         end
%         legend(reps, batchLabel);
%     end
    title([plotInfo.plotTitle ' cumulative distr'], 'Interpreter', 'tex', 'FontSize', 12);

end