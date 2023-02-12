% Driver script to generate specular and diffuse point clouds from
% measurements

savefile = 'geometry_test';

%%
x_dim = 200; y_dim = 200;
%x_dim = 190; y_dim = 140; % number of pixels along x and y dimension
bin_width = 8E-12; % bin time resolution in seconds
t0 = 14.2623e-9;% 14.2428e-9; % time corresponding to t=0

%% parameters for spot detector
win_size = 3; % 5 was for object scans, 3 was for big mirror/window
thresh_spot = 7; % 14 was for object scabs, 7 was for big mirror/window
plotDetectedSpots = true;

%% Glass processing parameters
refl_cal = 2.2316e4;

%%
data_dir = 'C:\Users\dolor\Remote-Flash\Experiments\Mirror Imaging\Paper\Github\Data\big_mirror\';
%data_dir = 'C:\Users\dolor\Remote-Flash\Experiments\window_disambiguation_test_031022\';

s = .257; % Laser-camera baseline separation
L = [s; 0; 0];

% Load laser pointing directions
lscanfile = [data_dir 'laser_scan_10x10.mat'];
%lscanfile = [data_dir 'laser_scan_14x10.mat'];
load(lscanfile, 'UV')

% Compute laser scanning azimuths and elevations
lthetas = acot(UV(:,1) ./ sqrt(1 + (UV(:,2)).^2));
lthetas(lthetas<0) = lthetas(lthetas<0)+pi;
lphis = atan(UV(:,2));

% Load detector scan angles
%cscanfile = [data_dir 'detector_scan_190x140.mat'];
cscanfile = [data_dir 'detector_scan_200x200.mat'];
load(cscanfile, 'thetaMap', 'phiMap')
thetaMap = flipud(reshape(thetaMap, x_dim, y_dim)');
phiMap = flipud(reshape(phiMap, x_dim, y_dim)');

% Load low level parameters extracted from raw data
params_file = [data_dir 'params.mat'];
load(params_file, 'pkbins', 'muvars', 'Evals', 'S')

%% Process frames and compute geometry

num_frames = size(UV, 1);
detections = struct('spots', {}, 'windows', {}, 'win_weights', {}, ...
    'spotcounts', {}, 'tofs', {}, 'thetas', {}, 'phis', {}, ...
    'ltheta', {}, 'lphi', {});

geometry = struct('D', {}, 'S1', {}, 'S2', {}, 'diffuse_first', {});
geometry(1).D = struct('pos', {}, 'r', {}, 'theta', {}, 'phi', {}, ...
   'tof', {}, 'ix', {}, 'rl', {}, 'counts', {}, 'reflectance', {}, ...
   'behind_window', {});
geometry(1).S1 = struct('pos', {}, 'n', {}, 'rl', {});
geometry(1).S2 = struct('pos', {}, 'n', {}, 'r', {}, 'thetas', {}, ...
   'phis', {}, 'tofs', {}, 'counts', {}, 'reflectance', {});

for ii = 1:num_frames
    
    disp(['Frame' num2str(ii) '.'])
    
    E_img = flipud(reshape(Evals(ii, :), x_dim, y_dim)');
    S_img = flipud(reshape(S(ii, :), x_dim, y_dim)');
    pkbin_img = flipud(reshape(pkbins(ii, :), x_dim, y_dim)');
    pkvar_img = flipud(reshape(muvars(ii, :), x_dim, y_dim)');

    frame_detections = processFrame(E_img, S_img, pkbin_img, pkvar_img, ...
        lthetas(ii), lphis(ii), thetaMap, phiMap, bin_width, t0, ...
        win_size, thresh_spot, plotDetectedSpots);
    
    detections(ii) = frame_detections;
  
    [pD, pS1, pS2, diffuse_first] = computeGeometryMirrorDisambiguation(detections(ii), s, refl_cal);
    
    geometry.D(ii) = pD;
    geometry.S1(ii) = pS1;
    geometry.S2(ii) = pS2;
    geometry.diffuse_first(ii) = diffuse_first;
    
end

save([data_dir savefile '.mat'], 'geometry', 'detections', 'thresh_spot', 'phiMap', 'thetaMap', 'lphis', 'lthetas', 'lscanfile', 'cscanfile');

%% Plot point cloud

disp('Plotting point cloud.')

% xlims = [-.5 1]; 
% ylims = [-1 .5]; 
% zlims = [0 2.5];

% xlims = [-.2 .6]; 
% ylims = [-.5 0]; 
% zlims = [.6 1.2];

xlims = [-1 1]; 
ylims = [-1.5 1.5]; 
zlims = [1.5 3];

normals_scale = 20; % Set to 0 to turn normals plotting off

save_pc = false; 

specularPointCloudPlotter(geometry, xlims, ylims, zlims, normals_scale)
view(-165, -60)
title('Point Cloud')

if save_pc
    saveas(gcf, [data_dir savefile '_pc.fig'])
end
