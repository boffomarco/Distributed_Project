clear all;
clc;


C = imread("Density100.png");
%imshow(C);

BW = rgb2gray(C);
figure,imshow(BW);

[X,Y] = meshgrid(1:100,1:100)
surf(X,Y,BW)

%{
% Set the figure
%figure(5)

% Set 3 robots randomly on the 
MatrixDimension = 100;
RobotsNumber = 5;
X = rand(RobotsNumber,2) * MatrixDimension;

 
%[v,c]=voronoi(X(:,1),X(:,2));

crs = [0,0;0,MatrixDimension;MatrixDimension,MatrixDimension;MatrixDimension,0];
[v,c]=VoronoiBounded(X(:,1),X(:,2), crs);

%polyarea(vx,vy)
%{
% Assign labels to the points.
nump = size(X,1);
plabels = arrayfun(@(n) {sprintf('X%d', n)}, (1:nump)');
hold on
Hpl = text(X(:,1), X(:,2), plabels, 'FontWeight', ...
      'bold', 'HorizontalAlignment','center', ...
      'BackgroundColor', 'none');


% Add a query point, P, at (5, 5).
P = [5 5];
plot(P(1),P(2), '*r');
text(P(1), P(2), 'P', 'FontWeight', 'bold', ...
     'HorizontalAlignment','center', ...
     'BackgroundColor', 'none');
hold off

%}
C = imread("Density100.png");
%imshow(C);

BW = rgb2gray(C);
figure,imshow(BW)


for i = 1:numel(c) %calculate the centroid of each cell
    [cx,cy] = PolyCentroid(BW,v(c{i},1),v(c{i},2));
    cx = min(MatrixDimension,max(0, cx));
    cy = min(MatrixDimension,max(0, cy));
    if ~isnan(cx) && inpolygon(cx,cy,uint8(crs(:,1)),uint8(crs(:,2)))
        X(:,1) = cx;  %don't update if goal is outside the polygon
        Y(:,1) = cy;
    end
end


%{
xq = zeros(1,1000);                                           % Random X-Coordinates
yq = zeros(1,1000);                                           % Random Y-Coordinates
for i=1:1000
    xq(i) = idivide(int16(i),100);
    yq(i) = rem( i , 100 );
end
xv = [1, 1, 50, 50, 1];                             % Polygon X-Coordinates
yv = [1, 50, 50, 1, 1];                             % Polygon Y-Coordinates

cumX = 0;
cumY = 0;
count = 0;
tot = 0;
for i=1:100

    for j=1:100
        

        [in,on] = inpolygon(i,j, xv,yv);                          % Logical Matrix
        inon = in | on;                                          % Combine ‘in’ And ‘on’
        if(inon)
            count = count + 1;
            tot = tot + uint32(BW(i,j));
            cumX = cumX + (j) * uint32(BW(i,j));
            cumY = cumY + (i) * uint32(BW(i,j));
        end
    end
end
X = (cumX/tot)
Y= (cumY/tot) 
%}


figure(1)
hold on                                       % Plot All Points
plot(cx, cy, 'gp')                                  % Overplot ‘inon’ Points
hold off





function [V,C]=VoronoiBounded(x,y, crs)
% VORONOIBOUNDED computes the Voronoi cells about the points (x,y) inside
% the bounding box (a polygon) crs.  If crs is not supplied, an
% axis-aligned box containing (x,y) is used.
bnd=[min(x) max(x) min(y) max(y)]; %data bounds
if nargin < 3
    crs=double([bnd(1) bnd(4);bnd(2) bnd(4);bnd(2) bnd(3);bnd(1) bnd(3);bnd(1) bnd(4)]);
end
rgx = max(crs(:,1))-min(crs(:,1));
rgy = max(crs(:,2))-min(crs(:,2));
rg = max(rgx,rgy);
midx = (max(crs(:,1))+min(crs(:,1)))/2;
midy = (max(crs(:,2))+min(crs(:,2)))/2;
% add 4 additional edges
xA = [x; midx + [0;0;-5*rg;+5*rg]];
yA = [y; midy + [-5*rg;+5*rg;0;0]];
[vi,ci]=voronoin([xA,yA]);
% remove the last 4 cells
C = ci(1:end-4);
V = vi;
% use Polybool to crop the cells
%Polybool for restriction of polygons to domain.
for ij=1:length(C)
        % thanks to http://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit
        % first convert the contour coordinate to clockwise order:
        [X2, Y2] = poly2cw(V(C{ij},1),V(C{ij},2));
        [xb, yb] = polybool('intersection',crs(:,1),crs(:,2),X2,Y2);
        ix=nan(1,length(xb));
        for il=1:length(xb)
            if any(V(:,1)==xb(il)) && any(V(:,2)==yb(il))
                ix1=find(V(:,1)==xb(il));
                ix2=find(V(:,2)==yb(il));
                for ib=1:length(ix1)
                    if any(ix1(ib)==ix2)
                        ix(il)=ix1(ib);
                    end
                end
                if isnan(ix(il))==1
                    lv=length(V);
                    V(lv+1,1)=xb(il);
                    V(lv+1,2)=yb(il);
                    ix(il)=lv+1;
                end
            else
                lv=length(V);
                V(lv+1,1)=xb(il);
                V(lv+1,2)=yb(il);
                ix(il)=lv+1;
            end
        end
        C{ij}=ix;
   
end
end


function [Cx,Cy] = PolyCentroid(BW, X,Y)
% POLYCENTROID returns the coordinates for the centroid of polygon with vertices X,Y
% The centroid of a non-self-intersecting closed polygon defined by n vertices (x0,y0), (x1,y1), ..., (xn?1,yn?1) is the point (Cx, Cy), where
% In these formulas, the vertices are assumed to be numbered in order of their occurrence along the polygon's perimeter, and the vertex ( xn, yn ) is assumed to be the same as ( x0, y0 ). Note that if the points are numbered in clockwise order the area A, computed as above, will have a negative sign; but the centroid coordinates will be correct even in this case.http://en.wikipedia.org/wiki/Centroid
% A = polyarea(X,Y)
cumX = 0;
cumY = 0;
count = 0;
tot = 0;
for i=1:1.0:100

    for j=1:1.0:100
        

        [in,on] = inpolygon(i,j, X,Y);                          % Logical Matrix
        inon = in | on;                                          % Combine ‘in’ And ‘on’
        if(inon)
            count = count + 1;
            tot = tot + uint8(BW(i,j));
            cumX = cumX + (j) * uint8(BW(i,j));
            cumY = cumY + (i) * uint8(BW(i,j));
        end
    end
end
Cx = (cumX/tot)
Cy = (cumY/tot)
end

%}