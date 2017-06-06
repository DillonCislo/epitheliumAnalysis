function out = makeDist(table, varargin)
    % makeDist(table, nBins)
    % makeDist(table, bins)
    %
    % function to make a distribution out of a table
    
    % bin number
    if nargin == 2 && numel(varargin{1}) == 1
        nBins = varargin{1};
    else
        nBins = 30;
    end
    
    % exclude NaN and Inf entries in table
    goodIdx = ~isnan(table) & abs(table) ~= Inf;
    newtable = table(goodIdx);
    
    % bins
    if nargin == 2 && numel(varargin{1}) > 1
        xbin = varargin{1};
    else
        avg = mean(newtable);
        sig = std(newtable);
        med = median(abs(newtable));
        
        % if sigma is huge because of outliers, use multiple of median
        % median of absolute value so distributions around zero dont go
        maxVal = min(avg + 3*sig, 15*med);
        xbin = 0:maxVal/(nBins-1):maxVal;
    end

    disp('FYI: makeDist excludes outliers with some hardcoded parameters');
    
    dist = hist(table(goodIdx), xbin);
    dist = dist/sum(dist(:));
    
    out = struct('dist', dist, 'bins', xbin, 'table', table);
end