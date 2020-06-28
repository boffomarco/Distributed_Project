clear all;
clc;

% initialize random generator
%rng('default')

%% Add the utility's functions folder
addpath('Utils')

%% VARIABLES
%declare global variabless
global abs_value_x; % N array, being the absolute x position of each robot
global abs_value_y; % N array, being the absolute y position of each robot
global counter;  % N array, being the counter for each robot position update
global mat_fig; % Figure to be covered by the robots (matrix form)

%define
N = 9; %number of robots
dim = 100; % Dimension of the picture square to cover
crs = [0,0;0,dim;dim,dim;dim,0]; % Polygon area definition
R = 25 + rand(N,1)*10; % Radius of the robots

%unCertainity = rand(N,1); % UnCertainty of the sensors on the robots
Covered = zeros(N,1); % Array used to store the Covered Area of each robot
TotCovered = []; % List used for the visualization of the total sum of Covered Area
StdCovered = []; % List used for the visualization of the std of Covered Area
iterations = 10000; % Number of executions of the algorithm
iter_tmp = 0; % store number of step to reach convergence for each image
% May include the name of the file inside images folder (equal to end_image to stop there)
image_i = 0; % Incremental to change image
end_image = 1; % Set number for final image +1 (>1 to test "moving images")
step = 0.7; % Maximum distance that a robot can move in a direction, 0.7 so to have maximum length of 1

abs_value_x = dim * rand(N,1)% Initial values of the robots
abs_value_y = dim * rand(N,1)% Initial values of the robots
counter = zeros(N,1);
pubs_abs_list = zeros(N,1);

% Debug variables
showPlot = true;
debug = false;


%% ROS NETWORK SETUP
% Initialize global node
%-------------------------------------
if not(ros.internal.Global.isNodeActive)     
    rosinit
end

% Topics and message types
%---------------------------------------
%source: http://docs.ros.org/melodic/api/geometry_msgs/html/msg/PointStamped.html
abs_value_topic = '/abs_value';
abs_value_msgs_type = 'geometry_msgs/PointStamped';
%source: http://docs.ros.org/api/nav_msgs/html/msg/OccupancyGrid.html
figure_topic = '/figure_image';
figure_msgs_type = 'nav_msgs/OccupancyGrid';

%List of robot nodes to be registered
%--------------------------------------
unreg_nodelist = cell(1,N);
unreg_tags = cell(1,N);
for i = 1:length(unreg_nodelist)
    node_name = char(strcat('robot_node_',string(i))); %name of node variable
    node_tag = char(strcat('/robot_',string(i))); %name of ros node (tag)
    unreg_tags(1,i) = {node_tag};
    unreg_nodelist(1,i) = {node_name};
end

%Initialize figure and multiple robot nodes, publishers and subscribers
%-------------------------------------------------------------------------
register = false; %reset flag
reg_nodelist = rosnode("list"); %list of registered nodes
if(length(reg_nodelist)==1) %check if no node has been register(apart from global node)
    %Register figure node
    figure_node = ros.Node('/figure');
    %Create publisher
    pub_figure = ros.Publisher(figure_node,figure_topic,figure_msgs_type);
    for i = 1:length(unreg_nodelist)       
        %Register nodes
        eval(sprintf('robot_node_%d = ros.Node(unreg_tags{i})', i));
        %Create publishers
        eval(sprintf('pub_abs_%d = ros.Publisher(robot_node_%d,abs_value_topic,abs_value_msgs_type)', i, i));
        %Save publishers in array
        eval(sprintf('pub_abs_list(%d) = pub_abs_%d',i,i));
        %Create subscribers
        eval(sprintf('sub_abs_%d = ros.Subscriber(robot_node_%d,abs_value_topic,abs_value_msgs_type,@newAbsCallback)', i, i));
        eval(sprintf('sub_fig_%d = ros.Subscriber(robot_node_%d,figure_topic,figure_msgs_type,@newFigCallback)', i, i));
    end    
else    
    for i = 1:length(unreg_nodelist)
        for j = 1:length(reg_nodelist)
            %Check if node alredy registered
            if (strcmp(unreg_tags{i}, reg_nodelist{j})) 
                register = true;
                break
            end
        end
        if not(register) %node not registered 
            %Register nodes
            eval(sprintf('robot_node_%d = ros.Node(unreg_tags{i})', i));
            %Create publishers
            eval(sprintf('pub_abs_%d = ros.Publisher(robot_node_%d,abs_value_topic,abs_value_msgs_type)', i, i));
            %Save publishers in array
            eval(sprintf('pub_abs_list(%d) = pub_abs_%d',i,i));
            %Create subscribers
            eval(sprintf('sub_abs_%d = ros.Subscriber(robot_node_%d,abs_value_topic,abs_value_msgs_type,@newAbsCallback)', i, i));
            eval(sprintf('sub_fig_%d = ros.Subscriber(robot_node_%d,figure_topic,figure_msgs_type,@newFigCallback)', i, i));
        end 
        register = false; %reset flag
    end
end

% Print ROS network
rosnode("list")


%% Visualization of initial points and Voronoi

close all
format compact

xrange = dim;
yrange = dim;

if showPlot
    verCellHandle = zeros(N,1);
    cellColors = cool(N);
    for i = 1:numel(abs_value_x) % color according to
        verCellHandle(i)  = patch(abs_value_x(i),abs_value_y(i),cellColors(i,:)); % use color i  -- no robot assigned yet
        hold on
    end
    pathHandle = zeros(N,1);    
    %numHandle = zeros(N,1);    
    for i = 1:numel(abs_value_x) % color according to
        pathHandle(i)  = plot(abs_value_x(i),abs_value_y(i),'-','color',cellColors(i,:)*.8);
    %    numHandle(i) = text(Px(i),Py(i),num2str(i));
    end
    goalHandle = plot(abs_value_x,abs_value_y,'+','linewidth',2);
    currHandle = plot(abs_value_x,abs_value_y,'o','linewidth',2);
    titleHandle = title(['o = Robots, + = Goals, Iteration ', num2str(0)]);
end
% End Visualization

%% Iterations
% For loop to iterate the algorithm
for iter = 1:iterations
    %[v,c]=VoronoiLimit(abs_value_x,abs_value_y, crs, false);
    [v,c]=VoronoiBounded(abs_value_x,abs_value_y, crs);
    
    % Print Delaunay triangulation
    if(debug && N >= 3)
        % Calculate the Delaunay triangulation
        t = delaunay ( abs_value_x, abs_value_y );
        % Display the Delaunay triangulation
        figure(2), clf, hold on;
        triplot ( t, abs_value_x, abs_value_y );
        title_string = sprintf ( 'Delaunay, step %d', iter );
        title ( title_string );
        axis equal
        %view ( 2 )
    end
    
    if showPlot 
        set(currHandle,'XData',abs_value_x,'YData',abs_value_y);%plot current position
        for i = 1:numel(abs_value_x) % color according to
            xD = [get(pathHandle(i),'XData'),abs_value_x(i)];
            yD = [get(pathHandle(i),'YData'),abs_value_y(i)];
            set(pathHandle(i),'XData',xD,'YData',yD, 'linewidth',2);%plot path position
     %       set(numHandle(i),'Position',[ abs_value_x(i),abs_value_y(i)]);
        end 
    end
    
    %% Update the figure
    % Do it just once per iteration instead of one for each robot at each iteration, since every robot will compute the same density function
    % [To be moved inside the iteration for each robot if we use multi agent simulation in MatLab]
    % [Is it possible that the robot build the density map by merging the information inside their Radius?]
    %{
    C = imread(strcat('Density',int2str(dim),'.png'));
    %imshow(C);
    BW = flipud(rgb2gray(C))';
    %figure,imshow(rgb2gray(C));
    if(debug)
        [X,Y] = meshgrid(1:dim,1:dim);
        figure(9)
        title_string = sprintf ( '3D figure at step %d', iter );
        title ( title_string );
        surf(X,Y,BW);
    end
    %}
    %--------------------------------------
    % Load image
    %img = imread(strcat('Density',int2str(dim),'.png'));
    img = imread(strcat('images/',int2str(image_i),'.png')); 
    mat_img = flipud(rgb2gray(img))';
    %imshow(mat_img)
    %Convert to occupancy map
    map = createOccupancyGrid(mat_img);
    %show(map)
    %Write msg
    msg_fig = rosmessage(figure_msgs_type);
    writeOccupancyGrid(msg_fig,map)
    %Send figure 
    send(pub_figure,msg_fig) % Sent from figure node
    pause(0.5) % Wait for message to update


    %% Iteration for each robot
    Rel = struct();
    for i = 1:N
        if(rand() < 2)% Randomly update position of robots
            [Cx,Cy, Rel, Covered(i)] = CoreAlgorithm(iter, mat_fig, i, abs_value_x, abs_value_y, R(i), dim, Rel);
            % Normalize movement
            if( double(int64(Cx)-abs_value_x(i)) > step)
                movX = step;
            else
                if( double(int64(Cx)-abs_value_x(i)) < - step)
                    movX = - step;
                else
                    movX = double(int64(Cx)-abs_value_x(i));
                    %movX = 0; % To avoid too many oscillations
                end
            end
            if( double(int64(Cy)-abs_value_y(i)) > step)
                movY = step;
            else
                if( double(int64(Cy)-abs_value_y(i)) < - step)
                    movY = - step;
                else
                    movY = double(int64(Cy)-abs_value_y(i));
                    %movY = 0; % To avoid too many oscillations
                end
            end
            
            %don't update if goal is outside the polygon
            if inpolygon(uint64(abs_value_x(i)+movX),uint64(abs_value_y(i)+movY),uint64(v(c{i},1)),uint64(v(c{i},2)))
                sendAbsValue(pub_abs_list, abs_value_msgs_type, abs_value_x(i) + movX, abs_value_y(i) + movY, i);
            end
        end
    end
    % End Iterations on the Robots
    
    %% Visualize Fitness values    
    tempTotCovered = 0;
    for i = 1:N 
        tempTotCovered = tempTotCovered + Covered(i); % Summation
    end
    TotCovered = [TotCovered tempTotCovered]; % Sum the coverage of each cell
    StdCovered = [StdCovered  std(Covered)]; % STD of the coverage of every cell
    figure(8), clf, hold on;
    title_string = sprintf ( 'Coverage & STD at step %d', iter );
    title ( title_string );
    yyaxis left
    ylabel('Total Coverage')
    plot(1:iter, TotCovered);
    yyaxis right
    ylabel('STD Coverage')
    plot(1:iter, StdCovered);
    
    %% Print Robots Network
    Xi = [];
    Yi = [];
    figure(3), clf, hold on;
    title_string = sprintf ( 'Robot positions, range & network at step %d', iter );
    title ( title_string );
    image('CData',mat_fig','XData',[0 dim],'YData',[0 dim])
    hold on
    for i=1:N
        plot(abs_value_x(i),abs_value_y(i),'wo','linewidth',2);
        %viscircles([abs_value_x(i) abs_value_y(i)], R(i), 'LineWidth',0.05);
    end
    for i=1:numel(Rel)
        if isfield(Rel(i), 'nei')
            plot(abs_value_x(i),abs_value_y(i), 'go','linewidth',2);
            viscircles([abs_value_x(i) abs_value_y(i)], R(i),'Color','white','LineStyle',':','linewidth',0.1);
            for j=1:numel(Rel(i).nei)
                Xi = [Xi abs_value_x(i) abs_value_x(Rel(i).nei(j)) NaN];
                Yi = [Yi abs_value_y(i)  abs_value_y(Rel(i).nei(j)) NaN];
            end
        end
    end
    plot(Xi',Yi', 'k', 'linewidth', 1);
    
    %% Show Voronoi Visualization
    if showPlot
        for i = 1:numel(c) % update Voronoi cells
            set(verCellHandle(i), 'XData',v(c{i},1),'YData',v(c{i},2));
        end
        set(titleHandle,'string',['o = Robots, + = Goals, Iteration ', num2str(iter,'%3d')]);
        set(goalHandle,'XData',abs_value_x,'YData',abs_value_y);%plot goal position
        
        axis equal
        axis([0,xrange,0,yrange]);
        drawnow
        %{
         if mod(iter,1) == 0
             %pause
             pause(0.001)
         end
        %}
    end
    
    %% End of computation in case of more than 100 steps or static image (count 2 near values in case of oscillation) (leave 2 steps to be sure)
    if((iter - iter_tmp > 100) || (numel(TotCovered)> 5 && ((TotCovered(end-4) == tempTotCovered && std(Covered) == StdCovered(end-4)) || (TotCovered(end-3) == tempTotCovered && std(Covered) == StdCovered(end-3)))))
        disp("Reached converged for image")
        image_i
        iter - iter_tmp
        iter_tmp = iter;
        image_i = image_i + 1;
        if(image_i == end_image)
            break;
        end
    end
end
% End Iterations

%% Visualize final Density function in 3D
figure(5)
title ( 'Final Density Function in 3D' );
[X,Y] = meshgrid(1:dim, 1:dim);
surf(X, Y, flipud(rgb2gray(img)));

%% Finish
% To check values 
abs_value_x;
abs_value_y;
counter;


