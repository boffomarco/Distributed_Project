function newFigCallback(~,msg) 
    global mat_fig;
    map = readOccupancyGrid(msg);
    mat_fig = occupancyMatrix(map); %Normalize values [0,1]   
end