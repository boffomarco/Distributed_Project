clear all;
clc;

addpath('Utils')

% Number of robots
num = 5;
% Dimension of area to cover
dim = 100;
% Radius of the robots
R = 30 + rand(num,1)*10;
% UnCertainty of the sensors on the robot
unCertainity = rand(num,1);
% Number of executions of the algorith
iterations = 1000;
% Algorithm to 
lloydsAlgorithm(R,num, dim, dim*rand(num,1),dim*rand(num,1), [0,0;0,dim;dim,dim;dim,0], iterations, true,false)



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

% Value used to get the covered area
Covered = zeros(n,1);
% Iteratively Apply LLYOD's Algorithm
for counter = 1:numIterations
    %[v,c]=VoronoiLimit(Px,Py, crs, false);
    [v,c]=VoronoiBounded(Px,Py, crs);
    
    
    if(debug && n >= 3)
        % Calculate the Delaunay triangulation
        t = delaunay ( Px, Py );
        % Display the Delaunay triangulation
        figure(2), clf, hold on;
        triplot ( t, Px, Py );
        title_string = sprintf ( 'Delaunay, step %d', counter );
        title ( title_string );
        axis equal
        %view ( 2 )
    end
    
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
    C = imread(strcat('Density',int2str(dim),'__.png'));
    %imshow(C);
    BW = flipud(rgb2gray(C))';
    %figure,imshow(rgb2gray(C));
    if(debug)
        [X,Y] = meshgrid(1:dim,1:dim);
        figure(3)
        surf(X,Y,BW);
    end
    
    Rel = struct();
    for i = 1:n %calculate the centroid of each cell
        if(rand()<2)% Randomly update position of robots
            [Cx,Cy, Rel, Covered(i)] = CoreAlgorithm(BW,i,Px,Py,R(i),dim, Rel);
            % Normalize movement
            if( double(int64(Cx)-Px(i)) > 0.75)
                movX = 0.75;
            else
                if( double(int64(Cx)-Px(i)) < -0.75)
                    movX = -0.75;
                else
                    movX = double(int64(Cx)-Px(i));
                end
            end
            if( double(int64(Cy)-Py(i)) > 0.75)
                movY = 0.75;
            else
                if( double(int64(Cy)-Py(i)) < -0.75)
                    movY = -0.75;
                else
                    movY = double(int64(Cy)-Py(i));
                end
            end
            
            %{
            if( double(int64(Cx)-Px(i)) > 0)
                movX = sqrt(double(int64(Cx)-Px(i)));
            else
                movX = -sqrt(norm(double(int64(Cx)-Px(i))));
            end
            if( double(int64(Cy)-Py(i)) > 0)
                movY = sqrt(double(int64(Cy)-Py(i)));
            else
                movY = -sqrt(norm(double(int64(Cy)-Py(i))));
            end
            %}
            
            %movX = double(int64(Cx)-Px(i));
            %movY = double(int64(Cy)-Py(i));
            %don't update if goal is outside the polygon
            if inpolygon(uint64(Px(i)+movX),uint64(Py(i)+movY),uint64(v(c{i},1)),uint64(v(c{i},2)))
                Px(i) = Px(i) + movX;
                Py(i) = Py(i) + movY;
            end
            %Px(i) = Cx;  
            %Py(i) = Cy;
        end
    end
    
    % Get the value of covered area
    TotCovered = 0;
    for i = 1:n %calculate the centroid of each cell
        TotCovered = TotCovered + Covered(i);
    end
    TotCovered
    
    Xi = [];
    Yi = [];
    figure(3), clf, hold on;
    image('CData',BW','XData',[0 dim],'YData',[0 dim])
    hold on
    for i=1:n
        plot(Px(i),Py(i),'w*')
        %viscircles([Px(i) Py(i)],R(i),'LineWidth',0.05);
    end
    for i=1:numel(Rel)
        if isfield(Rel(i),'nei')
            plot(Px(i),Py(i),'g*')
            viscircles([Px(i) Py(i)],R(i),'LineWidth',0.05);
            for j=1:numel(Rel(i).nei)
                Xi = [Xi Px(i) Px(Rel(i).nei(j)) NaN];
                Yi = [Yi Py(i)  Py(Rel(i).nei(j)) NaN];
            end
        end
    end
    plot(Xi',Yi');
    
    
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
figure(5)
[X,Y] = meshgrid(1:dim,1:dim);
surf(X,Y,flipud(rgb2gray(C)));
end