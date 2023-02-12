function detections = processFrame(E_img, S_img, pkbin_img, pkvar_img, ...
    ltheta, lphi, thetaMap, phiMap, bin_width, t0, ...
    win_size, thresh_spot, plotDetectedSpots)
% processFrame.m
% Connor Henley
% 1/31/2022
%
% Reads in frame information, detects visibile specularities and output
% relevant properties, including spot energy, time-of-flight, and angle of
% arrival.

% Inputs:
%   E_img = 2D image of detected photon counts - dim(y, x)
%   S_img = 2D image that specifies the pixels in which a backscattered 
%       pulse was detected (1 = detection, 0 = no detection)  - dim(y, x)
%   pkbin_img   = Image of fitted peak bins - dim(y, x) (x, y = spatial coordinates)
%   pkvar_img   = Image of peak bin fit variances - dim(y, x)
%   ltheta = Laser pointing azimuth (radians) - scalar
%   lphi = Laser point elevation (radians) - scalar
%   bin_width   = length of timing bin in seconds - scalar
%   t_0         = Time (in seconds) in receive window associated with zero 
%       time-of-flight - scalar
%   win_size   = window size for determining centroid of beam spot 
%       (sizes 3 and 5 are supported) - scalar
%   thresh_spot = "spottiness" threshold (see supplemental material) -
%       scalar
%   plotDetectedSpots = if true plots window and location of detected 
%       centroid - logical
%
% Outputs:
% detections = Structure that contains information about laser spots
%   identified in the frame.

% Store the laser pointing direction information 
detections.ltheta = ltheta;
detections.lphi = lphi;

% Detect spots in frame and compute time-of-flight associated with each
% spot
[spots, windows, win_weights, spotcounts] = spotDetector(E_img, S_img, win_size, thresh_spot, plotDetectedSpots);
tofs = binDetector(pkbin_img, pkvar_img, windows, bin_width, t0);

detections.spots = spots;
detections.windows = windows;
detections.win_weights = win_weights;
detections.spotcounts = spotcounts;
detections.tofs = tofs;

if ~isempty(tofs)
    % Angles of incidence for detected spots
    thetas = interp2(thetaMap, spots(:, 1), spots(:, 2));
    phis = interp2(phiMap, spots(:, 1), spots(:, 2));
    
    detections.thetas = thetas';
    detections.phis = phis';
else
    detections.thetas = [];
    detections.phis = [];
end

end

