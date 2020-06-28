%% Core Algorithm for hybrid approach
function [Cx,Cy, Rel, CoId] = CoreAlgorithm(iter, BW,index,Px,Py,R,dim,Rel)
% Set the borders to compute the Voronoi
minX = max(1,Px(index)-R);
maxX = min(dim,Px(index)+R);
minY = max(1,Py(index)-R);
maxY = min(dim,Py(index)+R);

% Store the original Density function to estimate the Covered Area
BWOriginal = BW;

% Visualization utils
Xi = [];
Yi = [];

% Get the neighbors, their position, their distance
VorX = [Px(index)];
VorY = [Py(index)];
Neighbours = [];
D = [];
for z=1:numel(Px)
    % calculate the norm of the distance
    dist = (Px(z)-Px(index))^2 + (Py(z)-Py(index))^2;
    % If the agent is inside the sensing area
    if  dist < R^2  && z ~= index
        % Save distance and position and index
        D = [D, sqrt(dist)];
        VorX = [VorX; Px(z)];
        VorY = [VorY; Py(z)];
        Neighbours = [Neighbours, z];
        
        % Show the density function affected by the neighbours of Robot 1
        if(index == 1)
            Xi = [Xi; Px(index); Px(z); NaN ];
            Yi = [Yi; Py(index); Py(z); NaN ];
        end
    end
end
% Save it to use later
Rel(index).nei = Neighbours;


% Update the density function for each detecter neighbour
for z=1:numel(Neighbours)
    % Get distance and compute euclidean norm to scale values
    [columnsInImage rowsInImage] = meshgrid(1:dim, 1:dim);
    circlePixels = (rowsInImage -  Px(Neighbours(z))).^2 + (columnsInImage - Py(Neighbours(z))).^2 <= R.^2;
    distanceImage = bwdist(~circlePixels);
    % Normalize it so that it's 1 in the center.
    distanceImage = distanceImage / max(distanceImage(:)); 
    %BW = BW - double(distanceImage*(255/numel(Neighbours))*((R-D(z))/R)); 
    % Add value to points which are close to neighbours which are far away
    BW = BW + double(distanceImage*(255/numel(Neighbours))*((D(z))/((D(z)+R)))); 
end
BW = max(BW,0.0); % Remove negative values

% Show the density function affected by the neighbours for Robot 1
if(index == 1)
    figure(4), clf, hold on;
    image('CData', BW', 'XData', [0 dim], 'YData', [0 dim])
    plot(Px(index),Py(index),'go','linewidth',2);
    title ( 'Robot(1): position, range, network & updated density function' );
    viscircles([Px(index) Py(index)],R,'LineWidth',0.05);
    plot(Xi,Yi, 'linewidth', 2);
end

% Calculate the Voronoi partition taking into account the neighbours
crs = [minX,minY;minX,maxY;maxX,maxY;maxX,minY];
[v,c] = VoronoiBounded(VorX,VorY, crs);

% Decreasing exponential when the number of Neighbours increases
expn = 2 * ( sqrt(numel(Neighbours)+1)/(numel(Neighbours)/1.555+0.2));
%expn = 2 * ( sqrt(numel(Neighbours)+1)/(numel(Neighbours)/1.555+0.2)) * (1+R/dim);

% Initialise the parameter to calculate the covered density area
CoId = 0;
% Initialise the parameters to compute the density function inside the
% voronoi partition
cumXV = 0;
totXV = 0;
cumYV = 0;
totYV = 0;
% Initialise the parameters to compute the density function inside the
% whole circle
cumX = 0;
totX = 0;
cumY = 0;
totY = 0;
%Iterate on subare of interest
for x=minX:1.0:maxX
    for y=minY:1.0:maxY
        % Iterate over values inside sensing area
        if (sqrt((Px(index)-x)^2 + (Py(index)-y)^2) <= R)
            % Compute check of Voronoi partitioning
            [in,on] = inpolygon(x,y,v(c{1},1),v(c{1},2)); % Logical Matrix
            inon = in | on;                               % Combine ‘in’ And ‘on’
            if(inon) % if point is inside Voronoi partitioning
                totXV = totXV + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
                totYV = totYV + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
                cumXV = cumXV + (x)  * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
                cumYV = cumYV + (y)  * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            end
            % also compute the points of the whole sensing area
            totX = totX + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            totY = totY + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            cumX = cumX + (x) * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            cumY = cumY + (y) * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            
            % Compute coverage of valued points inside sensing area
            CoId = CoId + double(int16(BWOriginal(uint64(x),uint64(y))));
        end
    end
end
% New centroids using only Voronoi partitioning (Scaled to have double weight wrt sensing)
CxV = (cumXV/(totXV))*3/2;
CyV = (cumYV/(totYV))*3/2;

% New centroids using whole sensing area (Scaled to have half weight wrt Voronoi)
Cx = (cumX/(totX))*1/2;
Cy = (cumY/(totY))*1/2;

% Update new centroids as hybrid approach
Cx = (CxV+Cx)/2;
Cy = (CyV+Cy)/2; 

end

