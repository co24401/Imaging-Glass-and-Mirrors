function naivePointCloudPlotter(geometry, xlims, ylims, zlims)
% naivePointCloudPlotter.m
% Connor Henley
% 3/5/2022
%
% Plot point cloud generated using naive_mirror_geometry_script

figure; 
for ii = 1:size(geometry.D, 2)
    hold on
    if ~isempty(geometry.D(ii).pos)
        inBeam = geometry.D(ii).inBeam;
        scatter3(geometry.D(ii).pos(1, inBeam), ...
            geometry.D(ii).pos(2, inBeam), ...
            geometry.D(ii).pos(3, inBeam), 10, 'b', 'filled')
        scatter3(geometry.D(ii).pos(1, ~inBeam), ...
            geometry.D(ii).pos(2, ~inBeam), ...
            geometry.D(ii).pos(3, ~inBeam), 10, 'b', '*')
    end
end

xlabel('X (cm)')
ylabel('Y (cm)')
zlabel('Z (cm)')

axis equal
grid on

xlim(xlims); 
ylim(ylims); 
zlim(zlims);

end

