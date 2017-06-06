function g = radialDistributionFunction(X,Y,xrange,yrange,rbins)
    % Calculate the radial distribution function
    %
    % g = radialDistributionFunction(X,Y,xrange,yrange,bins)
    %
    % X, Y:             Nx2 positions of particles
    % xrange, yrange:   range of particle position (x,y are X1, X2, don't confuse with X,Y)
    % bins:             radial bins (rmin:dr:rmax)
    % 
    % Definition (Chaikin & Lubensky p37):
    %
    % The pair distribution for a homogeneous system is:
    % g(x) =  <sum_{i ~= 0} delta(x -(x_i-x_0))>/<n>
    % Integrating this over a small volume gives the number of particles in
    % that volume, divided by the total density
    % at some x we measure the particles i for which x_i -x_0 = x, so which
    % are at distance x from the reference particle 0
    % When we integrate g over all all area, we should get the total number
    % of particles divided by the total density, so the total area.
    %
    % So we basically bin by:
    %
    % g(x) =  (A/N) (\int_a dA<sum_{i ~= 0} delta(x -(x_i-x_0))>)/a
    % where a is the area of the bin
    % 
    % We don't have ensemble averages, so we replace the average by
    % averaging over the origin:
    %
    % g(x) =  (A/N)(1/N)sum_j (\int_a dA sum_{i ~= j} delta(x -(x_i-x_j)))/a
    %
    % To bin radially and get the radial distribution function, we make a
    % an annulus between r and r+dr.
    %
    % To do correlation between distinct sets, we do
    % g(x) =  (A/Ny)(1/Nx) sum_j (\int_a dA sum_i delta(x -(x_i-y_j)))/a
    
    % area and range
    xmin = xrange(1);
    xmax = xrange(2);
    
    ymin = yrange(1);
    ymax = yrange(2);
    
    A = (xmax-xmin)*(ymax-ymin);
    
    % number of 'Y' particles
    Ny = size(Y,1);

    % radial bin related
    rmax = max(rbins);
    binArea = pi*(rbins(2:end).^2 - rbins(1:end-1).^2);
    
    % is this correlation or cross-correlation?
    symmetric = all(size(X)==size(Y)) && all(X(:)==Y(:));

    % keep track of cells in average
    Nx = 0;
    g = zeros([numel(rbins)-1 1]);
    
    for i = 1:size(X,1)

        % only average over cells that are far enough from image edge
        % i.e., only include this X in that case
        if ( X(i,1) > ymin + rmax && X(i,1) < ymax-rmax )...
            && ( X(i,2) > ymin + rmax && X(i,2) < xmax-rmax )

            dX = Y;
            dX(:,1) = Y(:,1) - X(i,1);
            dX(:,2) = Y(:,2) - X(i,2);

            if symmetric 
                dX(i,:) = [];
            end

            % density of Y cells as function of distance from X
            % throw out last bin because it only counts values matching the edge
            % (instead of an interval)
            P = histc(sqrt(sum(dX.^2,2)), rbins);
            P = P(1:end-1)./binArea';

            g = g + P;     
            Nx = Nx + 1;
        end
    end
    g = g*A/Ny;
    g = g/Nx;
end