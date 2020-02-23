function newAbsCallback(~,msg)    
    global abs_value; %array of Nx2 elements (being N #robots)
    abs_value(msg.Z,1) = msg.X; % Z = robot_id
    abs_value(msg.Z,2) = msg.Y;    
end
