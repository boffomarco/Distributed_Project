function map = createOccupancyGrid(img)
    fig_size = size(img);
    rows = fig_size(1);
    cols = fig_size(2);
    %Normalize values
    img_norm = mat2gray(img); %Normalize from 0-1
    img_norm_flip = flipud(img_norm); %Flip row order (due to the change of coordinates btw occupancy map and image mat)

    %Update occupancy grid 
    map = occupancyMap(rows,cols,1,'grid');    
    x = [1:cols]';
    y = ones([1,rows])';
    for this_row = img_norm_flip  
        updateOccupancy(map,[y x],this_row) 
        y = y + 1;
    end
    %show(map)
end