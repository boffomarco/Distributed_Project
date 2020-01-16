clear all;
clc;



C = imread("Density100.png");
imshow(C);

BW = rgb2gray(C);
figure,imshow(BW)

num = 3;
dim = 100;
lloydsAlgorithm(BW,num, dim, dim*rand(num,1),dim*rand(num,1), [0,0;0,dim;dim,dim;dim,0], 200, true)



function [Cx,Cy] = PolyCentroid(BW, X,Y)
% POLYCENTROID returns the coordinates for the centroid of polygon with vertices X,Y
% The centroid of a non-self-intersecting closed polygon defined by n vertices (x0,y0), (x1,y1), ..., (xn?1,yn?1) is the point (Cx, Cy), where
% In these formulas, the vertices are assumed to be numbered in order of their occurrence along the polygon's perimeter, and the vertex ( xn, yn ) is assumed to be the same as ( x0, y0 ). Note that if the points are numbered in clockwise order the area A, computed as above, will have a negative sign; but the centroid coordinates will be correct even in this case.http://en.wikipedia.org/wiki/Centroid
% A = polyarea(X,Y)
cumX = 0;
cumY = 0;
count = 0;
tot = 0;

%{
xq = zeros(1,10000);                                           % Random X-Coordinates
yq = zeros(1,10000);                                           % Random Y-Coordinates
for i=1:10000
    xq(i) = idivide(int16(i),100);
    yq(i) = rem( i , 100 );
end


[in,on] = inpolygon(xq,yq,  X, Y);                          % Logical Matrix
inon = in | on;        

idx = find(inon(:));                                        % Linear Indices Of �inon� Points
xcoord = xq(idx);                                           % X-Coordinates Of �inon� Points
ycoord = yq(idx);                                           % Y-Coordinates Of �inon� Points
figure(1)
plot(xq, yq, 'bp')                                          % Plot All Points
hold on
plot(X, Y, '-r')                                          % Plot Polygon
plot(xcoord, ycoord, 'gp')                                  % Overplot �inon� Points
hold off
%}
for i=1:1.0:100

    for j=1:1.0:100
        

        [in,on] = inpolygon(i,j, X,Y);                          % Logical Matrix
        inon = in | on;                                          % Combine �in� And �on�
        if(inon)
            count = count + 1;
            tot = tot + uint64(BW(i,j)+1);
            cumX = cumX + (i) * uint64(BW(i,j) + 1);
            cumY = cumY + (j) * uint64(BW(i,j) + 1);
        end
    end
end
Cx = (cumX/(tot))-1; 
Cy = (cumY/(tot))-1;

%{
Xa = [X(2:end);X(1)];
Ya = [Y(2:end);Y(1)];
A = 1/2*sum(X.*Ya-Xa.*Y); %signed area of the polygon

Cx = (1/(6*A)*sum((X + Xa).*(X.*Ya-Xa.*Y)));
Cy = (1/(6*A)*sum((Y + Ya).*(X.*Ya-Xa.*Y)));
%}

end


function [V,C]=VoronoiBounded(x,y, crs)
% VORONOIBOUNDED computes the Voronoi cells about the points (x,y) inside
% the bounding box (a polygon) crs.  If crs is not supplied, an
% axis-aligned box containing (x,y) is used.
bnd=[min(x) max(x) min(y) max(y)]; %data bounds
%{
if nargin < 3
    crs=double([bnd(1) bnd(4);bnd(2) bnd(4);bnd(2) bnd(3);bnd(1) bnd(3);bnd(1) bnd(4)]);
end
%}
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



function [Px, Py] = lloydsAlgorithm(BW,num,dim,Px,Py, crs, numIterations, showPlot)
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
%{
if nargin < 1   % demo mode
    showPlot = true;
    numIterations  = 200;
    xrange = 10;  %region size
    yrange = 5;
    n = 50; %number of robots  (changing the number of robots is interesting)
% Generate and Place  n stationary robots
    Px = 0.01*mod(1:n,ceil(sqrt(n)))'*xrange; %start the robots in a small grid
    Py = 0.01*floor((1:n)/sqrt(n))'*yrange;
    
%     Px = 0.1*rand(n,1)*xrange; % place n  robots randomly
%     Py = 0.1*rand(n,1)*yrange;
    
    crs = [ 0, 0;    
        0, yrange;
        1/3*xrange, yrange;  % a world with a narrow passage
        1/3*xrange, 1/4*yrange;
        2/3*xrange, 1/4*yrange;
        2/3*xrange, yrange;
        xrange, yrange;
        xrange, 0];
    
    for i = 1:numel(Px)  
        while ~inpolygon(Px(i),Py(i),crs(:,1),crs(:,2))% ensure robots are inside the boundary
            Px(i) = rand(1,1)*xrange; 
            Py(i) = rand(1,1)*yrange;
        end
    end
else
    xrange = max(crs(:,1));
    yrange = max(crs(:,2));
    n = numel(Px); %number of robots  
end
%}
n = num;

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
    
    if showPlot
        set(currHandle,'XData',Px,'YData',Py);%plot current position
        for i = 1:numel(Px) % color according to
            xD = [get(pathHandle(i),'XData'),Px(i)];
            yD = [get(pathHandle(i),'YData'),Py(i)];
            set(pathHandle(i),'XData',xD,'YData',yD);%plot path position
     %       set(numHandle(i),'Position',[ Px(i),Py(i)]);
        end 
    end
    
    for i = 1:numel(c) %calculate the centroid of each cell
        [cx,cy] = PolyCentroid(BW,v(c{i},1),v(c{i},2))
        cx = min(xrange,max(0, cx));
        cy = min(yrange,max(0, cy));
        if ~isnan(cx) && inpolygon(uint64(cx),uint64(cy),uint64(crs(:,1)),uint64(crs(:,2)))
            Px(i) = cx;  %don't update if goal is outside the polygon
            Py(i) = cy;
        end
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