function specularPointCloudPlotter(geometry, xlims, ylims, zlims, n_scale)
% specularPointCloudPlotter.m
% Connor Henley
% 3/11/2022
%
% Plot point cloud generated using MirrorGeometryScript

if n_scale ~= 0
    plot_normals = true;
else
    plot_normals = false;
end

figure; 
for ii = 1:size(geometry.D, 2)
    hold on
    if ~isempty(geometry.D(ii).pos)
        behind_window = geometry.D(ii).behind_window;
        if geometry.diffuse_first(ii)
            D = geometry.D(ii).pos(:, ~behind_window);
            scatter3(D(1, :),  D(2, :), D(3, :), 10, 'b', 'filled')
            Dout = geometry.D(ii).pos(:, behind_window);
            scatter3(Dout(1, :),  Dout(2, :), Dout(3, :), 10, 'b', 'o')
        else
            behind_window = geometry.D(ii).behind_window;
            D = geometry.D(ii).pos(:, ~behind_window);
            scatter3(D(1, :),  D(2, :), D(3, :), 10, 'b', '*')
            Dout = geometry.D(ii).pos(:, behind_window);
            scatter3(Dout(1, :),  Dout(2, :), Dout(3, :), 10, 'b', 'o')
        end
    end
    
    if ~isempty(geometry.S1(ii).pos)
        if geometry.diffuse_first(ii)
            scatter3(geometry.S1(ii).pos(1), ...
                geometry.S1(ii).pos(2), ...
                geometry.S1(ii).pos(3), 10, 'g', 'filled')
        else
            scatter3(geometry.S1(ii).pos(1), ...
                geometry.S1(ii).pos(2), ...
                geometry.S1(ii).pos(3), 10, 'g', '*')
        end
        if plot_normals
            quiver3(geometry.S1(ii).pos(1), ...
                geometry.S1(ii).pos(2), ...
                geometry.S1(ii).pos(3), ...
                geometry.S1(ii).n(1)/n_scale, ...
                geometry.S1(ii).n(2)/n_scale, ...
                geometry.S1(ii).n(3)/n_scale, 0, 'g')
        end
    end
    
    if ~isempty(geometry.S2(ii).pos)
        if geometry.diffuse_first(ii)
            scatter3(geometry.S2(ii).pos(1, :), ...
                geometry.S2(ii).pos(2, :), ...
                geometry.S2(ii).pos(3, :), 10, 'r', 'filled')
        else
            scatter3(geometry.S2(ii).pos(1,:), ...
                geometry.S2(ii).pos(2,:), ...
                geometry.S2(ii).pos(3,:), 10, 'r', '*')
        end
        if plot_normals
            quiver3(geometry.S2(ii).pos(1,:), ...
                geometry.S2(ii).pos(2,:), ...
                geometry.S2(ii).pos(3,:), ...
                geometry.S2(ii).n(1,:)/n_scale, ...
                geometry.S2(ii).n(2,:)/n_scale, ...
                geometry.S2(ii).n(3,:)/n_scale, 0, 'r')
        end
    end
end

xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

axis equal
grid on

xlim(xlims); 
ylim(ylims); 
zlim(zlims);

end

