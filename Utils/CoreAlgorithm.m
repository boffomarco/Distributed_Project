
function [Cx,Cy, Rel, Covered] = CoreAlgorithm(BW,index,Px,Py,R,dim,Rel)
% Set the borders to compute the Voronoi
minX = max(1,Px(index)-R);
maxX = min(dim,Px(index)+R);
minY = max(1,Py(index)-R);
maxY = min(dim,Py(index)+R);

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
        BW = BW - double(distanceImage*255*((R-sqrt(D))/R));
        
        VorX = [VorX; Px(z)];
        VorY = [VorY; Py(z)];
        Neighbours = [Neighbours, z];
        
        
        Xi = [Xi; Px(index); Px(z); NaN ];
        Yi = [Yi; Py(index); Py(z); NaN ];
    end
end
Rel(index).nei = Neighbours;

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

expn = (2 * ( numel(Neighbours) + 1)) ^2;
%expn = 2;

% Initialise the parameter to calculate the covered density area
Covered = 0;
% Initialise the parameters to compute the density function inside the
% voronoi partition
cumX = 0;
totX = 0;
cumY = 0;
totY = 0;
% Initialise the parameters to compute the density function inside the
% whole circle
cumXT = 0;
totXT = 0;
cumYT = 0;
totYT = 0;
for x=minX:1.0:maxX
    for y=minY:1.0:maxY
        if (sqrt((Px(index)-x)^2 + (Py(index)-y)^2) <= R)
            [in,on] = inpolygon(x,y,v(c{1},1),v(c{1},2));                          % Logical Matrix
            inon = in | on;                                          % Combine ‘in’ And ‘on’
            if(inon)
                totX = totX + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
                totY = totY + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
                cumX = cumX + (x)  * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
                cumY = cumY + (y)  * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            end
            
            totXT = totXT + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            totYT = totYT + (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            cumXT = cumXT + (x) * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            cumYT = cumYT + (y) * (double(int16(BW(uint64(x),uint64(y)))+1)/256)^(expn);
            
            Covered = Covered + double(int16(BW(uint64(x),uint64(y))));
        end
    end
end

Cx = (cumX/(totX));
Cy = (cumY/(totY));

Cxt = (cumXT/(totXT));
Cyt = (cumYT/(totYT));

Cx = (Cx+Cxt)/2;
Cy = (Cy+Cyt)/2;
end

