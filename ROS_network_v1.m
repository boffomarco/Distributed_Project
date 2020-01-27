%% VARIABLES
N = 10; %number of robots



%% ROS NETWORK SETUP
% Initialize global node
%-------------------------------------
if not(ros.internal.Global.isNodeActive)     
    rosinit
end
%clear all

% Topics and message types
%---------------------------------------
%source: http://docs.ros.org/api/geometry_msgs/html/msg/Point.html
abs_value_topic = '/abs_value';
abs_value_msgs_type = 'geometry_msgs/Point';
%source: http://docs.ros.org/api/nav_msgs/html/msg/OccupancyGrid.html
%http://docs.ros.org/api/sensor_msgs/html/msg/Image.html
figure_topic = '/figure_image';
figure_msgs_type = 'std_msgs/String';

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
        %Create subscribers
        eval(sprintf('sub_abs_%d = ros.Subscriber(robot_node_%d,abs_value_topic,abs_value_msgs_type)', i, i));
        eval(sprintf('sub_fig_%d = ros.Subscriber(robot_node_%d,figure_topic,figure_msgs_type)', i, i));
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
            %Create subscribers
            eval(sprintf('sub_abs_%d = ros.Subscriber(robot_node_%d,abs_value_topic,abs_value_msgs_type)', i, i));
            eval(sprintf('sub_fig_%d = ros.Subscriber(robot_node_%d,figure_topic,figure_msgs_type)', i, i));
        end 
        register = false; %reset flag
    end
end

% Print ROS network
rosnode("list")

%% SEND & RECEIVE MESSAGES (EXAMPLE)
% Update figure image
msg_fig = rosmessage(figure_msgs_type);
msg_fig.Data = 'Figure_image';

%Send figure to topic 
send(pub_figure,msg_fig) % Sent from figure node
pause(1) % Wait for message to update

%Retrieve figure
sub_fig_1.LatestMessage.Data
sub_fig_2.LatestMessage.Data
sub_fig_3.LatestMessage.Data

%Send absolute position value
msg_abs_val = rosmessage(abs_value_msgs_type);

msg_abs_val.X = 12;
msg_abs_val.Y = 21;
send(pub_abs_1,msg_abs_val)
pause(1) % Wait for message to update

%Retrieve abs position
% sub_abs_2.LatestMessage.Y
% sub_abs_3.LatestMessage.Z


