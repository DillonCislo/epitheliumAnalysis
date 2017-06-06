function [ClockSorted,Index] = ClockWiseSortCell(cell, vertices)

        NewCellLength = length(cell);
        
        CM = mean(vertices(cell,:),1);
        
        CMVertices = vertices(cell,:);
        CMVertices(:,2) = CMVertices(:,2) - CM(2);
        CMVertices(:,1) = CMVertices(:,1) - CM(1);
        
        [Sorted,Index] = sort(cart2pol(CMVertices(:,1),CMVertices(:,2))); 
        
        % Now we actually sort counter clock, just invert.
        Index = Index(length(Index):-1:1);
        ClockSorted = vertices(cell(Index),:); 
end