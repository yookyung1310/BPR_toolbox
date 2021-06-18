function mmconverted = unitconversion(au_obs, coeff,distance_prev,distance,base)
% 2018.03.19 KY updated
% this function converts pupil size (au) to pupil diameter (mm) with
% approxiamtation to circle
% See the explaining slide. 
% coeff        : computed coefficient based on the dots measured
% distnace_prev: the distance between the dot and eyetracker lens when the coefficient computed (cm)
% distnace     : present distance between the dot and eyetracker lens (cm)
% base         : present horizontal distance between eyetracker and eye (cm)


% First, adjust the measured orthogonal projection as parallel
theta1 = acos(base/distance);
% NOTE: gazeHeigth is relative height from the fixation level, so it should be negative if eye level is lower
% than fixation level
gazeHeight = 0; % 10 % from eye level to gaze level (fixation point level) (cm)
theta2 = atan(gazeHeight/distance);

theta = theta1+theta2;
au_real = au_obs/cos(theta);

% Second, adjust the coefficient based on present distance.
coeff_adj  = coeff*(distance/distance_prev); 

% now convert it!
mmconverted = coeff_adj*2*sqrt(au_real/pi);


end