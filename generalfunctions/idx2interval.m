function interval = idx2interval(idx)
    % convert list of integers into integer intervals
    % 
    % interval = idx2interval(idx)
    %
    % idx:      list of integers
    % interval: Nx2 list of intervals
    
    vals = 1:max(idx);
    mask = false([max(idx), 1]);
    mask(idx) = true;

    ends = mask - circshift(mask,-1) == 1;
    starts = mask - circshift(mask,+1) == 1;
    
    % first and last element are not necessarily picked up by above trick
    if min(idx) == 1
        starts(1) = 1;
    end
    ends(end) = 1;

    interval = [vals(starts)', vals(ends)']; 
end