function newFigCallback(~,msg) 
    global mat_fig;
    map = readOccupancyGrid(msg);
    mat_fig = occupancyMatrix(map) .* double(255.0); %Normalize values [0,255] 
end