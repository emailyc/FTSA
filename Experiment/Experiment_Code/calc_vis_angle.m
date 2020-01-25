% Xiaoqin (Feb 2015)
%
% To use function:
% screenwidth = horiozontal or vertical screen width (mm)
% resolution  = horizontal or vertical screen resolution (pixels)
% viewingdist = viewing distance (mm)
%
% Note: Dimensions of screenwidth and resolution have to correspond
%
% To convert pixels to degrees:
% degree = deg_per_pixel*no_of_pixels
% 
% To convert degrees to pixels:
% pixels = (1/deg_per_pixel)*degrees

function pixel_per_deg = calc_vis_angle(screenwidth, resolution, viewingdist)

pixel_per_deg = pi * (resolution) / atan(screenwidth/viewingdist/2) / 360;

end