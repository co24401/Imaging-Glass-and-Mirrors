% Driver script to implement specular surface mapping using multi-beam
% illumination.

noise_gate_start = 1;
noise_gate_stop = 1000;
min_pk_cts_threshold = 0.8;

x_dim = 200; y_dim = 200; % number of pixels along x and y dimension
bin_width = 8E-12; % bin time resolution in seconds
t0 = 14.2623e-9; % time corresponding to t=0

win_size = 3; %5;
thresh_spot = 7;
plotDetectedSpots = true;

%% Load data and relevant parameters

data_dir = 'C:\Users\dolor\Remote-Flash\Experiments\Mirror Imaging\Paper\Github\Data\window_multibeam\';
dataCubeFile = [data_dir 'dataCube_summed.mat'];
load(dataCubeFile,  'dataCube_summed')

replicafile = 'C:\Users\dolor\Remote-Flash\Experiments\Mirror Imaging\Paper\Github\Data\replica\replica.mat';
load(replicafile, 'mu_fit_LUT', 'pkbias_LUT')

lscanfile = [data_dir 'laser_scan_11x9.mat'];
load(lscanfile, 'UV')

lthetas = acot(UV(:,1) ./ sqrt(1 + (UV(:,2)).^2));
lthetas(lthetas<0) = lthetas(lthetas<0)+pi;
lphis = atan(UV(:,2));

lthetas = lthetas';
lphis = lphis';

% Load detector scan angles
cscanfile = [data_dir 'detector_scan_200x200.mat'];
load(cscanfile, 'thetaMap', 'phiMap', 'uMap', 'vMap', 'num_u', 'num_v', 'u_lims', 'v_lims')
thetaMap = flipud(reshape(thetaMap, x_dim, y_dim)');
phiMap = flipud(reshape(phiMap, x_dim, y_dim)');
uMap = flipud(reshape(uMap, x_dim, y_dim)');
vMap = flipud(reshape(vMap, x_dim, y_dim)');

c = 299792458; % speed of light (m/s)
s = .257; % Laser-camera baseline separation
L = [s; 0; 0]; % source position
in_beam_tolerance = .000075;% Maximum cosine distance of point from beam ray
out_beam_buffer = .000075;
%   for point to classified as lying along laser beam path
% .000075; worked for mirror dataset

%% Spot Extraction

% Extract per-pixel maps of detected photons, peak ToF bin, etc.
[pkbins, mubins, muvars, Evals, Evars, ...
    w2s, S] = datacubeStats(dataCube_summed, ...
    noise_gate_start, noise_gate_stop, min_pk_cts_threshold, ...
    mu_fit_LUT, pkbias_LUT);

E_img = flipud(reshape(Evals, x_dim, y_dim)');
S_img = flipud(reshape(S, x_dim, y_dim)');
pkbin_img = flipud(reshape(pkbins, x_dim, y_dim)');
pkvar_img = flipud(reshape(muvars, x_dim, y_dim)');

% Detect spots and extract times of flight.
[spots, windows, win_weights, spotcounts] = spotDetector(E_img, S_img, win_size, thresh_spot, plotDetectedSpots);
tofs = binDetector(pkbin_img, pkvar_img, windows, bin_width, t0);

% Compute angles of arrival
spotus = interp2(reshape(uMap, x_dim, y_dim), spots(:, 1), spots(:, 2))';
spotvs = interp2(reshape(vMap, x_dim, y_dim), spots(:, 1), spots(:, 2))';
thetas = acot(spotus ./ sqrt(1 + (spotvs).^2));
thetas(thetas<0) = thetas(thetas<0)+pi;
phis = atan(spotvs);

% Naive position estimates
r_xc = 0.5*(c^2 * tofs.^2 - s^2)./( c*tofs - s*cos(thetas) );
dir_points = [cos(thetas); sin(thetas).*sin(phis); sin(thetas).*cos(phis)];
X = r_xc.*dir_points; % Interesting to plot all of these.  You can see the distortion that occurs when you use the wrong source position.
dir_beams = [cos(lthetas); sin(lthetas).*sin(lphis); sin(lthetas).*cos(lphis)];

% Beam intersection
dir_XL_est = (X - L) ./ vecnorm(X - L, 2, 1);
sepBX = 1 - dir_beams'*dir_XL_est;

[mnxs, ixs] = min(sepBX, [], 2);  % Find spot that is closest to each beam
inds = ixs(mnxs < in_beam_tolerance); % Closest spot to each beam must be within angular distance threshold
spotsInBeam = false(size(tofs));
spotsInBeam(ixs) = true;

XiB = X(:, spotsInBeam); % Apparent positions of spots along beam vectors.
uiB = spotus(spotsInBeam);
viB = spotvs(spotsInBeam);

spotsOutBeam = ~spotsInBeam;
uoB = spotus(spotsOutBeam);
voB = spotvs(spotsOutBeam);

% Use measured Z position of in-beam points to interpolate Z position of
% out of beam points as a function of incidence direction.
Z = scatteredInterpolant(uiB',viB',XiB(3, :)'); % Default interp method is linear
zoB = Z(uoB', voB')';
XoB = zoB.*[uoB; voB; ones(size(uoB))];

r_xc_oB = vecnorm(XoB, 2, 1);

toB = tofs(spotsOutBeam);

uu = linspace(u_lims(2), u_lims(1), num_u);
vv = linspace(v_lims(1), v_lims(2), num_v);
figure; imagesc(uu, vv, reshape(Evals, num_u, num_v)')
set(gca, 'XDir', 'reverse')
set(gca, 'YDir', 'normal')
set(gca, 'ColorScale', 'log')
colorbar
caxis([5 3000])
hold on; scatter(uiB, viB, 'og')
hold on; scatter(uoB, voB, 'or')

%% Mirror source position localization

% Initialize mirror source position with true source position
X_0_init = [0; 0; 0];

t = 0.5; % Damping constant for Newton step (1 = undamped).
epsilon = 1e-10; %0.0000001; % Stopping criterion.
max_iter = 1000;

k = 1000; % Num RANSAC iters
n = 4; % Number of maybeInliers
d = 10; % Num. alsoInliers
min_g2 = 0.025;  % Inlier threshold
allInds = 1:size(XoB, 2);

X_0_best = X_0_init;
msqerr_best = Inf;
errIx_best = [];
msqerr_record = [];

for ii = 1:k

    maybeInInds = randsample(size(XoB, 2), n);
    maybeIn = ismember(allInds, maybeInInds);
    otherPts = ~ismember(allInds, maybeIn);

    [X_0_test,~] = ...
        mirrorSourceFit(XoB(:, maybeIn), toB(maybeIn), ...
        X_0_init, t, epsilon, max_iter, false);

    dx = X_0_test(1) - XoB(1,:);
    dy = X_0_test(2) - XoB(2,:);
    dz = X_0_test(3) - XoB(3,:);
    dt = c*toB - r_xc_oB;
    g = dx.^2 + dy.^2 + dz.^2 - dt.^2;

    alsoIn = g.^2 < min_g2 & otherPts;
    numAlsoIn = sum(alsoIn);
    disp([num2str(numAlsoIn) ' also-inliers.'])

    if numAlsoIn > d
        errIx = alsoIn | maybeIn;
        msqerr = mean(g(errIx).^2);

        if msqerr < msqerr_best
            msqerr_best = msqerr;
            msqerr_record(end+1) = msqerr;
            X_0_best = X_0_test;
            errIx_best = errIx;
        end
    end

end

[X_0,~] = ...
    mirrorSourceFit(XoB(:, errIx_best), toB(errIx_best), ...
    X_0_best, t, epsilon, max_iter, true);

%% Compute mirror plane

r_0 = norm(X_0);
n_mirror = -X_0/r_0;
d_mirror = -r_0/2;
P_mirror = [n_mirror; -d_mirror];

nd_dot = (n_mirror'*dir_beams);

dir_mirror_beams =  dir_beams - 2*nd_dot.*n_mirror;

figure; quiver3(zeros(1, size(dir_beams,2)), zeros(1, size(dir_beams,2)), zeros(1, size(dir_beams,2)), dir_beams(1,:), dir_beams(2, :), dir_beams(3,:), 'r')
hold on; quiver3(X_0(1)*ones(1, size(dir_beams,2)), X_0(2)*ones(1, size(dir_beams,2)), X_0(3)*ones(1, size(dir_beams,2)), dir_mirror_beams(1,:), dir_mirror_beams(2, :), dir_mirror_beams(3,:), 'b')
hold on; scatter3(X_0(1)/2, X_0(2)/2, X_0(3)/2, 80, 'kp')
legend({'True beams', 'Mirrored Beams', 'Mirror Point'})
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
axis equal
grid on
view(-165, -60)

%% Ray-plane intersections

nd_dot_mirror = (n_mirror'*dir_mirror_beams);

cos_deltas_0 = (X_0'*dir_points)/r_0;

r_xc_mirror = (1/2)*(r_0^2 - c^2*tofs.^2)./(r_0.*cos_deltas_0 - c*tofs);

X_mirror = r_xc_mirror.*dir_points;  % Apparent position of 2B points

dir_points_X_0 = (X_mirror - X_0)./vecnorm(X_mirror - X_0, 2, 1);

intersections_CM = dir_points.*d_mirror./(n_mirror'*dir_points);
intersections_0M = X_0 + dir_points_X_0.*...
    (d_mirror - n_mirror'*X_0)./ (n_mirror'*dir_points_X_0);

in_front_of_mirror = n_mirror'*X_mirror-d_mirror > 0;
behind_mirror = n_mirror'*X_mirror-d_mirror < 0;

mirror_points_C = intersections_CM(:, ~spotsInBeam & behind_mirror);
mirror_points_0 = intersections_0M(:, ~spotsInBeam & in_front_of_mirror);

figure; scatter3(X(1,spotsInBeam), X(2, spotsInBeam), X(3, spotsInBeam), 'filled');
hold on; scatter3(X_mirror(1,~spotsInBeam), X_mirror(2, ~spotsInBeam), X_mirror(3, ~spotsInBeam), 'filled')
hold on; scatter3(mirror_points_C(1,:), mirror_points_C(2,:), mirror_points_C(3,:), 'filled')
hold on; scatter3(mirror_points_0(1,:), mirror_points_0(2,:), mirror_points_0(3,:), 'filled')
axis equal
legend({'1B 3B spots', 'False 2B pos.', '2B spots, from camera' , '2B spots, from mirror source'})

%% Flip 3B points across mirror

mirror_points = [mirror_points_0 mirror_points_C];
uMirror = mirror_points(1,:)./mirror_points(3,:);
vMirror = mirror_points(2,:)./mirror_points(3,:);
k_hull = convhull(uMirror, vMirror);

behindMirror = inpolygon(uiB, viB, uMirror(k_hull), vMirror(k_hull));
X1B = XiB(:, ~behindMirror);
X3B = XiB(:, behindMirror);
X3B = X3B + 2*d_mirror*n_mirror - 2*(n_mirror'*X3B).*n_mirror;

%% Plot point cloud

figure; scatter3(X1B(1,:), X1B(2,:), X1B(3,:), 10, 'b', 'filled')
hold on; scatter3(X3B(1,:), X3B(2,:), X3B(3,:), 10, 'b', '*')
hold on; scatter3(mirror_points_C(1,:), mirror_points_C(2,:), mirror_points_C(3,:), 10, 'r', 'filled')
hold on; scatter3(mirror_points_0(1,:), mirror_points_0(2,:), mirror_points_0(3,:), 10, 'g', 'filled')
axis equal
legend({'1B', 'Flipped 3B', '2B: rays from C' , '2B: rays from X_0'})
xlim([-1 1]); ylim([-1.5 1.5]); zlim([0 3]);
view(-160, -58)
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
