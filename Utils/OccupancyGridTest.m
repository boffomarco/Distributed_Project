%% BUILD MAT FIGURE TEST
C = imread("Density100.png");
fig_size = size(C);
rows = fig_size(1);
cols = fig_size(2);
%Normalize values
C_norm = mat2gray(C); %Normalize from 0-1
C_norm_flip = flipud(C_norm); %Flip row order (due to the change of coordinates btw occupancy map and image mat)

%Update occupancy grid 
map = occupancyMap(rows,cols,1,'grid');
x = [1:cols]';
y = ones([1,rows])';
for this_row = C_norm_flip  
    updateOccupancy(map,[y x],this_row) 
    y = y + 1;
end
show(map)
%Convert back to 
%mat_figure = round(255*occupancyMatrix(map));
mat_figure = occupancyMatrix(map); %Normalize values [0,1]
%imshow(mat_figure)