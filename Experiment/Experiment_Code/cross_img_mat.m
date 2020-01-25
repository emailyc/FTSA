% Xiaoqin 11 Feb 2015
% Creates a image matrix of a fixation cross
% cross_size = size of the cross
% thickness = 


function fix_cross = cross_img_mat(cross_size, thickness, cross_col, backgrd_col)

fix_cross(1:cross_size, 1:cross_size, 1:3) = cross_col;
fix_cross(1:cross_size/2-round(thickness/2), 1:cross_size/2-round(thickness/2), :) = backgrd_col;
fix_cross(cross_size/2+round(thickness/2):end, 1:cross_size/2-round(thickness/2), :) = backgrd_col;
fix_cross(1:cross_size/2-round(thickness/2), cross_size/2+round(thickness/2):end, :) = backgrd_col;
fix_cross(cross_size/2+round(thickness/2):end, cross_size/2+round(thickness/2):end, :) = backgrd_col;
fix_cross = uint8(fix_cross);
