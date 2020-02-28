clear all
%%GLOBAL VARIABLES
global chat_ex;

%% MAIN

% Create a master node 
rosshutdown
rosinit
rostopic list

%Initialize multiple nodes.
node1 = ros.Node('/test_node_3');
node2 = ros.Node('/test_node_4');

%Use these nodes to perform separate operations and send separate messages
%A message published by node1 can be accessed by a subscriber running in node2.
pub = ros.Publisher(node1,'/chatter','std_msgs/String');
sub = ros.Subscriber(node2,'/chatter','std_msgs/String');

%Create sub with callback (callback function will be triggered everytime a
%new message is published
%Callbacks
sub_cb = rossubscriber('/chatter',@callbackExample);

msg = rosmessage('std_msgs/String');
msg.Data = 'Message from Node 1';

%Send a message from node1. The subscriber attached to node2 will receive the message
send(pub,msg) % Sent from node 1
pause(1) % Wait for message to update
sub.LatestMessage

%Callback variable update
fprintf('This is the value received by callback: %s\n\n', chat_ex);

%Clear the ROS network of publisher, subscriber, and nodes. 
% Delete the Core object to shut down the ROS master.
rosshutdown

%% CALLBACK FUNCTION

function callbackExample(~,msg)
    global chat_ex;
    chat_ex = msg.Data;    
end

%% USEFUL FUNCTIONS

% %Show details of a ros message
% ros_msg_details = rosmessage('std_msgs/String');

% %Show details of a message
% details = showdetails(msg);




