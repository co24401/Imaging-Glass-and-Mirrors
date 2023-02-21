% Driver script that generates point cloud from measurements assuming that 
% all detections correspond to one-bounce light transport paths.  

data_dir = 'C:\Users\dolor\Remote-Flash\Experiments\window_disambiguation_test_031022\'; %'C:\Users\dolor\Remote-Flash\Experiments\window_collection_peepers_01822\';
detections_file = 'geometry_mda_strictinbeam.mat'; %'test_window_tolerant_geometry.mat';
savefile = 'naive_geometry';

load([data_dir detections_file], 'detections')

s = .257; % Laser-camera baseline separation

geometry(1).D = struct('pos', {}, 'r', {}, 'theta', {}, 'phi', {}, ...
    'tof', {}, 'counts', {}, 'inBeam', {});


for ii = 1:length(detections)

    pD = computeGeometryNaively(detections(ii), s);
    geometry.D(ii) = pD;
    
end

save([data_dir savefile '.mat'], 'geometry', 'detections');
%% Plot point cloud

disp('Plotting point cloud.')

xlims = [-1 1.25]; 
ylims = [-1.5 .5]; 
zlims = [0 3.25];

save_pc = true;

naivePointCloudPlotter(geometry, xlims, ylims, zlims)
view(-165, -60)
title('Point Cloud')

if save_pc
    saveas(gcf, [data_dir savefile '_pc.fig'])
end