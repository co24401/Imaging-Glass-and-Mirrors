% Spot position figure in ECCV paper

geom_file = 'C:\Users\dolor\Remote-Flash\Experiments\window_collection_nopeepers_011722\geometry_february.mat';

load(geom_file, 'geometry')

diffuse_first = logical(geometry.diffuse_first);
D_df = [geometry.D(diffuse_first).pos];
D_sf = [geometry.D(~diffuse_first).pos];
S2_df = [geometry.S2(diffuse_first).pos];
S2_sf = [geometry.S2(~diffuse_first).pos];

figure; imagesc(uu, vv, reshape(Evals, num_u, num_v)')
set(gca, 'XDir', 'reverse')
set(gca, 'YDir', 'normal')
set(gca, 'ColorScale', 'log')
colorbar
caxis([5 3000])
hold on; scatter([D_df(1,:)./D_df(3,:) S2_sf(1,:)./S2_sf(3,:)], [D_df(2,:)./D_df(3,:) S2_sf(2,:)./S2_sf(3,:)], 'og')
hold on; scatter([D_sf(1,:)./D_sf(3,:) S2_df(1,:)./S2_df(3,:)], [D_sf(2,:)./D_sf(3,:) S2_df(2,:)./S2_df(3,:)], 'or')
set(gca, 'FontSize', 14)