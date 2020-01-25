function [fix_cross, rect_cross] = CreateFixationCross(window, param)

% Cross size
cross_size_pix = round(param.ppd*param.cross_size);
if mod(cross_size_pix, 2) ~= 0
    cross_size_pix = cross_size_pix + 1;
end

% Cross thickness
cross_thick_pix = round(param.ppd*param.cross_thickness);
if mod(cross_thick_pix, 2) ~= 0
    cross_thick_pix = cross_thick_pix + 1;
end

fix_cross_temp = CreateCrossMat(cross_size_pix, cross_thick_pix, param.stim_col(1), param.cross_bg_col(1));

% Define rectangle which we want the cross to be in
rect_cross = [param.center_x-floor(cross_size_pix/2) param.center_y-floor(cross_size_pix/2) param.center_x+floor(cross_size_pix/2) param.center_y+floor(cross_size_pix/2)];

% Make fixation cross
fix_cross = Screen('MakeTexture', window, fix_cross_temp);