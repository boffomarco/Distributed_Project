function newAbsCallback(~,msg)    
    global abs_value_x; %array of Nx2 elements (being N #robots)
    global abs_value_y;
    global counter;
    abs_value_x(msg.Point.Z) = msg.Point.X; % Z = robot_id
    abs_value_y(msg.Point.Z) = msg.Point.Y;   
    counter(msg.Point.Z) = msg.Header.Seq; % Seq number, it increments automatically
end
