function [V,C] = randomVoronoiLattice2(a, noise, mask)
    % generate a 2D random lattice as Voronoi cell of a noisy dual
    %
    % [V,C] = randomVoronoiLattice2(a, noise, mask)
    %
    % [V,C]:    Voronoi lattice:
    % V:        Nx2 vertices
    % C:        cell array of vertices
    %
    % a:        cell size (~lattice constant)
    % noise:    magnitude of noise added to dual lattice
    % mask:     bw image
    %
    % Idse Heemskerk, 2015
    
    ySize = size(mask,1);
    xSize = size(mask,2);

    [X,Y] = meshgrid(-a:a:xSize+a, -a:a*sqrt(3):ySize+a);
    X = cat(1, X(:), X(:) + a/2);
    Y = cat(1, Y(:), Y(:) + a*sqrt(3)/2);

    X = X + noise*(rand(size(X)) - 0.5);
    Y = Y + noise*(rand(size(X)) - 0.5);

    % first crop at mask image borders
    ind1 = X < xSize & Y < ySize & X > 0 & Y > 0;

    % then crop outside actual mask
    ind = sub2ind(size(mask), 1+floor(Y(ind1)), 1+floor(X(ind1)));
    ind1(ind1) = mask(ind);
    
    % make voronoi
    [vvor,cvor] = voronoin([X,Y]);
    cvor = cvor(ind1); 

    % trim cells
    cc = {};
    for i = 1:length(cvor)
        % v(1) == Inf, if the cell has no vertex at infinity add it to cc
        % if all vertices are also within the approximate grid boundaries
        if all(cvor{i}~=1)
            cc = [cc cvor(i)];
        end
    end

    % remove unused vertices
    vidx = false(length(vvor),1);
    % cc{i} contain indices of vertices
    for i = 1:length(cc)
        vidx(cc{i}) = 1;
    end
    vertices = vvor(vidx,1:2);
    vshift = zeros(length(vvor),1);
    vshift(vidx) = 1:sum(vidx);
    for i = 1:length(cc)
        cc{i} = vshift(cc{i})';
    end
    
    vertices = [vertices 0*vertices(:,1)];
    V = vertices;
    C = cc;
end