function [V,C] = randomVoronoiLattice(N, a, noise)
    % generate a 2D random lattice as Voronoi cell of a noisy dual
    %
    % [V,C]:    Voronoi lattice:
    % V:        Nx2 vertices
    % C:        cell array of vertices
    %
    % N:        controls number of cells (but indirectly)
    % a:        anisotropy, as ratio of long over short axis
    % noise:    magnitude of noise added to dual lattice
    %
    % lattice spacing 1
    %
    % Idse Heemskerk, 2015

    if nargin == 0
        N = 3; 
        a = 1; % anisotropy
        noise = 0.7;
    end

    % generate a triangular grid
    [X,Y] = meshgrid(-2*N:2*N, -N:N);
    Y = Y*a*sqrt(3); 
    X = cat(1, X(:), X(:) + 1/2);
    Y = cat(1, Y(:), Y(:) + a*sqrt(3)/2);

    % approximate grid boundaries
    Lx = 2*N - 1;
    Ly = 2*N - 1;

    % scatter(X(:),Y(:),'.r');
    % axis equal;

    % add some noise
    X = X + noise*(rand(size(X)) - 0.5);
    Y = Y + noise*(rand(size(X)) - 0.5);

    % the voronoi tesselation with tension points as seeds is a vertex model
    % that  admits the delaunay tension net
    [vvor,cvor] = voronoin([X,Y]);

%     voronoi(X,Y, 'r');
%     axis equal;

    % trim cells
    cc = {};
    for i = 1:length(cvor)
        % v(1) == Inf, if the cell has no vertex at infinity add it to cc
        % if all vertices are also within the approximate grid boundaries
        if all(cvor{i}~=1)...
                && all(vvor(cvor{i},1) < Lx) && all(vvor(cvor{i},1) > -Lx)...
                && all(vvor(cvor{i},2) > -Ly) && all(vvor(cvor{i},2) < Ly)

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
    
    V = vertices;
    C = cc;
end