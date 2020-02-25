function sendAbsValue(pub_id_list, msg_type, abs_x,abs_y,robot_id)
    %Create message container
    msg_abs_val = rosmessage(msg_type);
    %Fill the messge
    msg_abs_val.Point.Z = robot_id; %robot id
    msg_abs_val.Point.X = abs_x;
    msg_abs_val.Point.Y = abs_y;
    %Send the message
    send(pub_id_list(robot_id),msg_abs_val)
    pause(1) % Wait for message to update
end