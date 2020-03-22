
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

% Get the neighbors, their position, and update the density function
VorX = [Px(index)];
VorY = [Py(index)];
Neighbours = [];
for z=1:numel(Px)
    D = (Px(z)-Px(index))^2 + (Py(z)-Py(index))^2;
    if  D < R^2  && z ~= index
        [columnsInImage rowsInImage] = meshgrid(1:dim, 1:dim);
        circlePixels = (rowsInImage -  Px(z)).^2 + (columnsInImage - Py(z)).^2 <= R.^2;
        distanceImage = bwdist(~circlePixels);
        % Normalize it so that it's 1 in the center.
        distanceImage = distanceImage / max(distanceImage(:));
        %BW = BW - double(distanceImage*255*((R-sqrt(D))/R));
        BW = BW - double(distanceImage*(R-sqrt(D)));
        
        VorX = [VorX; Px(z)];
        VorY = [VorY; Py(z)];
        Neighbours = [Neighbours, z];
        
        % Show the density function affected by the neighbours
        if(index == 1)
            Xi = [Xi; Px(index); Px(z); NaN ];
            Yi = [Yi; Py(index); Py(z); NaN ];
        end
    end
end
Rel(index).nei = Neighbours;
BW = max(BW,0.0); % Remove negative values

% Calculate the Voronoi partition taking into account the neighbours
crs = [minX,minY;minX,maxY;maxX,maxY;maxX,minY];
[v,c] = VoronoiBounded(VorX,VorY, crs);

% Show the density function affected by the neighbours
if(index == 1)
    figure(4), clf, hold on;
    image('CData', BW', 'XData', [0 dim], 'YData', [0 dim])
    plot(Px(index),Py(index),'g*')
    title ( 'Robot(1): position, range, network & updated density function' );
    viscircles([Px(index) Py(index)],R,'LineWidth',0.05);
    plot(Xi,Yi);
end

% Decreasing exponential when the number of Neighbours increases
expn = 2 * ( sqrt(numel(Neighbours)+1)/(numel(Neighbours)+1) + 1);
% Start the first 50 iterations to just cover homogeneusly the area, 
% then cover the density image
if(iter < 50)
    expn = 1 / expn;
end

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
for x=minX:1.0:maxX
    for y=minY:1.0:maxY
        if (sqrt((Px(index)-x)^2 + (Py(index)-y)^2) <= R)
            [in,on] = inpolygon(x,y,v(c{1},1),v(c{1},2));                          % Logical Matrix
            inon = in | on;                                          % Combine ‘in’ And ‘on’
            if(inon)
                
                totXV = totXV + (double(int16(BWOriginal(uint64(x),uint64(y)))+1)/256)^(expn);
                totYV = totYV + (double(int16(BWOriginal(uint64(x),uint64(y)))+1)/256)^(expn);
                cumXV = cumXV + (x)  * (double(int16(BWOriginal(uint64(x),uint64(y)))+1)/256)^(expn);
                cumYV = cumYV + (y)  * (double(int16(BWOriginal(uint64(x),uint64(y)))+1)/256)^(expn);
                
            end
            
            totX = totX + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            totY = totY + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            cumX = cumX + (x) * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            cumY = cumY + (y) * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            
            CoId = CoId + double(int16(BWOriginal(uint64(x),uint64(y))));
        end
    end
end

CxV = (cumXV/(totXV));
CyV = (cumYV/(totYV));

Cx = (cumX/(totX));
Cy = (cumY/(totY));

Cx = (CxV+Cx)/2;
Cy = (CyV+Cy)/2; 

end

