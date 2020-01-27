% Create a master node
rosinit
rostopic list

%Initialize multiple nodes.
node1 = ros.Node('/test_node_3');
node2 = ros.Node('/test_node_4');

%Use these nodes to perform separate operations and send separate messages
%A message published by node1 can be accessed by a subscriber running in node2.
pub = ros.Publisher(node1,'/chatter','std_msgs/String');
sub = ros.Subscriber(node2,'/chatter','std_msgs/String');

msg = rosmessage('std_msgs/String');
msg.Data = 'Message from Node 1';

%Send a message from node1. The subscriber attached to node2 will receive the message
send(pub,msg) % Sent from node 1
pause(1) % Wait for message to update
sub.LatestMessage

%Clear the ROS network of publisher, subscriber, and nodes. 
% Delete the Core object to shut down the ROS master.
rosshutdown

%% USEFUL FUNCTIONS

% %Show details of a ros message
% ros_msg_details = rosmessage('std_msgs/String');

% %Show details of a message
% details = showdetails(msg);


