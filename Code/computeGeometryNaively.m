function pD = computeGeometryNaively(detections, s)
% computeGeometry2.m
% Connor Henley
% 3/5/2022
%
% Compute 3D points from detections.  Assumes all detected spots are 1B
% returns.  Records whether or not the return is in the beam, however,

% Inputs:
% detections = Structure that information about spots detected in
% the frame, such as time and angle of arrival for each spot.
% s = Baseline seperation between source (L) and camera aperture (origin).
%
% Outputs:
% pD = Structure that holds information about point on diffusely reflecting
% surface (D)

c = 299792458; % speed of light (m/s)
L = [s; 0; 0]; % source position
in_beam_tolerance = 0.01; %0.001; % Maximum cosine distance of point from beam ray
%   for point to classified as lying along laser beam path

% Load detection params to de-clutter code
tofs = detections.tofs;
thetas = detections.thetas;
phis = detections.phis;
spotcounts = detections.spotcounts;
ltheta = detections.ltheta;
lphi = detections.lphi;

pD = struct('pos', [], 'r', [], 'theta', [], 'phi', [], ...
    'tof', [], 'counts', [], 'inBeam', []);

if ~isempty(tofs)
    
    r_xc = 0.5*(c^2 * tofs.^2 - s^2)./( c*tofs - s*cos(thetas) );
    X = r_xc.*[cos(thetas); sin(thetas).*sin(phis); sin(thetas).*cos(phis)];
    dir_point_L = [cos(ltheta); sin(ltheta)*sin(lphi); sin(ltheta)*cos(lphi)];
    
    dir_XL_est = (X - L) ./ vecnorm(X - L, 2, 1);
    sepFromBeam = (1 - dir_point_L'*dir_XL_est);
    inBeam = sepFromBeam < in_beam_tolerance;
    
    pD.pos = X;
    pD.r = r_xc;
    pD.theta = thetas;
    pD.phi = phis;
    pD.tof = tofs; 
    pD.counts = spotcounts;
    pD.inBeam = inBeam;
    
end

end

