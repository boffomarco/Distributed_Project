clear all;
clc;

num = 20;
dim = 100;
R = 50;
lloydsAlgorithm(R,num, dim, dim*rand(num,1)/10+dim/2,dim*rand(num,1)/10+dim/2, [0,0;0,dim;dim,dim;dim,0], 200, true,false)

function [Cx,Cy] = CoreAlgorithm(BW,index,Px,Py,R,dim)

minX = min(0,Px(index)-R);
maxX = max(dim,Px(index)+R);
minY = min(0,Py(index)-R);
maxY = max(dim,Py(index)+R);

VorX = [Px(index)];
VorY = [Py(index)];
for z=1:numel(Px)
    if (((Px(z)-Px(index))^2 + (Py(z)-Py(index))^2)^(1/2) < R) && z ~= index
        VorX = [VorX; Px(z)];
        VorY = [VorY; Py(z)];
    end
end

crs = [minX,minY;minX,maxY;maxX,maxY;maxX,minY];
[v,c]=VoronoiBounded(VorX,VorY, crs);

cumX = 0;
totX = 0;
cumY = 0;
totY = 0;

cumXT = 0;
totXT = 0;
cumYT = 0;
totYT = 0;
for x=minX:1.0:maxX
    for y=minY:1.0:maxY
        if  uint64(x)>0 &&  uint64(x)<=100 &&  uint64(y)>0 &&  uint64(y)<=dim
            [in,on] = inpolygon(x,y,v(c{1},1),v(c{1},2));                          % Logical Matrix
            inon = in | on;                                          % Combine ‘in’ And ‘on’
            if(inon)
                if ((Px(index)-x)^2 + (Py(index)-y))^(1/2) < R
                    totX = totX + ((dim - ((Px(index)-x)^2)^(1/2)) )^2 * uint64(BW(uint64(x),uint64(y))+1)^2;
                    totY = totY + ((dim - ((Py(index)-y)^2)^(1/2)))^2 * uint64(BW(uint64(x),uint64(y))+1)^2;
                    cumX = cumX + (x) * ((dim - ((Px(index)-x)^2)^(1/2)) )^2 * (uint64(BW(uint64(x),uint64(y)) + 1))^2;
                    cumY = cumY + (y) * ((dim - ((Py(index)-y)^2)^(1/2)))^2 * (uint64(BW(uint64(x),uint64(y)) + 1))^2;
                end

                totXT = totXT + ((dim - ((Px(index)-x)^2)^(1/2)) )^2 * uint64(BW(uint64(x),uint64(y))+1)^2;
                totYT = totYT + ((dim - ((Py(index)-y)^2)^(1/2)))^2 * uint64(BW(uint64(x),uint64(y))+1)^2;
                cumXT = cumXT + (x) * ((dim - ((Px(index)-x)^2)^(1/2)) )^2 * (uint64(BW(uint64(x),uint64(y)) + 1))^2;
                cumYT = cumYT + (y) * ((dim - ((Py(index)-y)^2)^(1/2)))^2 * (uint64(BW(uint64(x),uint64(y)) + 1))^2;
            end
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

function [Cx,Cy] = PolyCentroid(BW,Xcrs,Ycrs, X,Y, dim)
% POLYCENTROID returns the coordinates for the centroid of polygon with vertices X,Y
% The centroid of a non-self-intersecting closed polygon defined by n vertices (x0,y0), (x1,y1), ..., (xn?1,yn?1) is the point (Cx, Cy), where
% In these formulas, the vertices are assumed to be numbered in order of their occurrence along the polygon's perimeter, and the vertex ( xn, yn ) is assumed to be the same as ( x0, y0 ). Note that if the points are numbered in clockwise order the area A, computed as above, will have a negative sign; but the centroid coordinates will be correct even in this case.http://en.wikipedia.org/wiki/Centroid
% A = polyarea(X,Y)
cumX1 = 0;
cumY1 = 0;
totX1 = 0;
totY1 = 0;
cumX2 = 0;
cumY2 = 0;
totX2 = 0;
totY2 = 0;
for i=1:1.0:dim
    for j=1:1.0:dim
        [in,on] = inpolygon(i,j, Xcrs,Ycrs);                          % Logical Matrix
        inon = in | on;                                          % Combine ‘in’ And ‘on’
        if(inon)
            totX1 = totX1 + ((dim - ((X-i)^2)^(1/2)) )^2 * uint64(BW(i,j)+1)^2;
            totY1 = totY1 + ((dim - ((Y-j)^2)^(1/2)))^2 * uint64(BW(i,j)+1)^2;
            cumX1 = cumX1 + (i) * ((dim - ((X-i)^2)^(1/2)) )^2 * (uint64(BW(i,j) + 1))^2;
            cumY1 = cumY1 + (j) * ((dim - ((Y-j)^2)^(1/2)))^2 * (uint64(BW(i,j) + 1))^2;
        end
        
        totX2 = totX2 + ((dim - ((X-i)^2)^(1/2)) )^2 * (uint64(BW(i,j)+1))^2;
        totY2 = totY2 + ((dim - ((Y-j)^2)^(1/2)) )^2  * (uint64(BW(i,j)+1))^2;
        cumX2 = cumX2 + (i) * ((dim - ((X-i)^2)^(1/2)))^2 * (uint64(BW(i,j)+1))^2;
        cumY2 = cumY2 + (j) * ((dim - ((Y-j)^2)^(1/2)))^2 * (uint64(BW(i,j)+1))^2;

    end
end
Cx1 = (cumX1/(totX1)); 
Cy1 = (cumY1/(totY1));

Cx2 = (cumX2/(totX2)); 
Cy2 = (cumY2/(totY2));


Cx = (Cx1+Cx2)/2; 
Cy = (Cy1+Cy2)/2; 

end


function [V,C]=VoronoiBounded(x,y, crs)
% VORONOIBOUNDED computes the Voronoi cells about the points (x,y) inside
% the bounding box (a polygon) crs. 
bnd=[min(x) max(x) min(y) max(y)]; %data bounds

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



function [Px, Py] = lloydsAlgorithm(R,n,dim,Px,Py, crs, numIterations, showPlot,debug)
% LLOYDSALGORITHM runs Lloyd's algorithm on the particles at xy positions 
% (Px,Py) within the boundary polygon crs for numIterations iterations
% showPlot = true will display the results graphically.  
% 
% Lloyd's algorithm starts with an initial distribution of samples or
% points and consists of repeatedly executing one relaxation step:
%   1.  The Voronoi diagram of all the points is computed.
%   2.  Each cell of the Voronoi diagram is integrated and the centroid is computed.
%   3.  Each point is then moved to the centroid of its Voronoi cell.
%
% Inspired by http://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit
% Requires the Polybool function of the mapping toolbox to run.
%
% Run with no input to see example.  To initialize a square with 50 robots 
% in left middle, run:
%lloydsAlgorithm(0.01*rand(50,1),zeros(50,1)+1/2, [0,0;0,1;1,1;1,0], 200, true)
%
% Made by: Aaron Becker, atbecker@uh.edu
% close all
format compact
% initialize random generator in repeatable fashion
sd = 20;
rng(sd)

xrange = dim;
yrange = dim;
%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if showPlot
    verCellHandle = zeros(n,1);
    cellColors = cool(n);
    for i = 1:numel(Px) % color according to
        verCellHandle(i)  = patch(Px(i),Py(i),cellColors(i,:)); % use color i  -- no robot assigned yet
        hold on
    end
    pathHandle = zeros(n,1);    
    %numHandle = zeros(n,1);    
    for i = 1:numel(Px) % color according to
        pathHandle(i)  = plot(Px(i),Py(i),'-','color',cellColors(i,:)*.8);
    %    numHandle(i) = text(Px(i),Py(i),num2str(i));
    end
    goalHandle = plot(Px,Py,'+','linewidth',2);
    currHandle = plot(Px,Py,'o','linewidth',2);
    titleHandle = title(['o = Robots, + = Goals, Iteration ', num2str(0)]);
end
%%%%%%%%%%%%%%%%%%%%%%%% END VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Iteratively Apply LLYOD's Algorithm
for counter = 1:numIterations
    %[v,c]=VoronoiLimit(Px,Py, crs, false);
    [v,c]=VoronoiBounded(Px,Py, crs);
    
    % Calculate the Delaunay triangulation
    t = delaunay ( Px, Py );
    % Display the Delaunay triangulation
    %{
    figure(2), clf, hold on;
    triplot ( t, Px, Py );
    title_string = sprintf ( 'Delaunay, step %d', counter );
    title ( title_string );
    axis equal
    view ( 2 )
    %}
    
    if showPlot 
        set(currHandle,'XData',Px,'YData',Py);%plot current position
        for i = 1:numel(Px) % color according to
            xD = [get(pathHandle(i),'XData'),Px(i)];
            yD = [get(pathHandle(i),'YData'),Py(i)];
            set(pathHandle(i),'XData',xD,'YData',yD);%plot path position
     %       set(numHandle(i),'Position',[ Px(i),Py(i)]);
        end 
    end
    
    
    % Update the figure    
    C = imread("Density100.png");
    %imshow(C);
    BW =rgb2gray(C);
    %figure,imshow(BW);
    if(debug)
        [X,Y] = meshgrid(1:100,1:100);
        figure(3)
        surf(X,Y,BW);
    end
    
    for i = 1:n %calculate the centroid of each cell
        [Cx,Cy] = CoreAlgorithm(BW,i,Px,Py,R,dim);
        movX = double(int64(Cx)-Px(i));
        movY = double(int64(Cy)-Py(i));
        if inpolygon(uint64(Px(i)+movX),uint64(Py(i)+movY),uint64(v(c{i},1)),uint64(v(c{i},2)))
            Px(i) = Px(i) + movX;
            Py(i) = Py(i) + movY;
        end
        %Px(i) = Cx;  %don't update if goal is outside the polygon
        %Py(i) = Cy;
    end
    
    if showPlot
        for i = 1:numel(c) % update Voronoi cells
            set(verCellHandle(i), 'XData',v(c{i},1),'YData',v(c{i},2));
        end
        set(titleHandle,'string',['o = Robots, + = Goals, Iteration ', num2str(counter,'%3d')]);
        set(goalHandle,'XData',Px,'YData',Py);%plot goal position
        
        axis equal
        axis([0,xrange,0,yrange]);
        drawnow
        %{
         if mod(counter,1) ==0
             %pause
             pause(0.001)
         end
        %}
    end
    
end

end