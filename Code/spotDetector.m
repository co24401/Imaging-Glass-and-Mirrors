%  Detects laser spots in data and extracts their angle of arrival and 
%  intensity (in photon counts).
%
% Inputs:
%   E_img = 2D image of detected photon counts - dim(y, x)
%   S_img = 2D image that specifies the pixels in which a backscattered 
%       pulse was detected (1 = detection, 0 = no detection)  - dim(y, x)
%   win_size   = window size for determining centroid of beam spot 
%       (sizes 3 and 5 are supported) - scalar
%   thresh_spot = "spottiness" threshold (see supplemental material) -
%       scalar
%   plotDetectedSpots = if true plots window and location of detected 
%       centroid - logical
%
% Output:
%   spots = array of indices where spot is detected - dim(s, 2) list (s = # of detected spots)
%   windows = Boundaries of spot detection window - dim(s, 4) list 
%   win_weights = Normalized per-pixel intensities within spot fitting
%       window - dim(s, 3, 3)
%   spotcounts = Total number of photon detections associated with detected
%       spot - dim(s)
%
%% NOTES
%   * make sure there is sufficient space at top and bottom of image (i.e.
%     no spots detected within few pixels of edge.
%   * If unsure choose following values:
%       win_size   = 3
%       thresh_pct = 7

function [spots, windows, win_weights, spotcounts] = spotDetector(E_img, S_img, win_size, thresh_spot, plotDetectedSpots)
%% make sure you have image processing toolbox
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
    error('You do not have image processing toolbox needed to run code');
end

x_dim = size(E_img, 2);
y_dim = size(E_img, 1);

%% convolve with laplacian to detect edges
% Need to update this to adapt to arbitrary window size.  Currently only
% produces 3x3 and 5x5 filters
if win_size == 3
    h_spot = fspecial('laplacian', 0.2);
    spot_img = -1*conv2(E_img, h_spot);
    spot_img = spot_img(2:end-1, 2:end-1); % reduce dimension from convolution
    h_avg = fspecial('average', 3);
    avg_img = conv2(E_img, h_avg);
    avg_img = avg_img(2:end-1, 2:end-1); % reduce dimension from convolution
else
    h_spot = [-1 -1 -1 -1 -1; -1 1 2 1 -1; -1 2 4 2 -1; -1 1 2 1 -1; -1 -1 -1 -1 -1];
    spot_img = conv2(E_img, h_spot);
    spot_img = spot_img(3:end-2, 3:end-2); % reduce dimension from convolution
    h_avg = fspecial('average', 5);
    avg_img = conv2(E_img, h_avg);
    avg_img = avg_img(3:end-2, 3:end-2); % reduce dimension from convolution
end



%% threshold to detect blobs/spots
binary_img = S_img & (spot_img./avg_img > thresh_spot);
labeled_img = bwlabel(binary_img, 8);
numSpots = max(max(labeled_img));

if numSpots == 0
    warning('No spots detected.');
end

%% extract spot windows, determine centroids
x_mins = []; x_maxs = [];
y_mins = []; y_maxs = [];
spots = [];
windows = zeros(numSpots, 4);
spotcounts = zeros(1, numSpots);
maxcorrvals = zeros(1, numSpots);
win_weights = zeros(win_size, win_size, numSpots);
for s=1:numSpots
    %% create mask for a given spot
    maskedImg = labeled_img==s;

    %% extract highest intensity pixel in spot
    maskedImg = maskedImg .* spot_img;
    max_val = max(max(maskedImg));
    maxcorrvals(s) = max_val;
    [y_cent, x_cent] = find(maskedImg == max_val);

    %% determine window range
    dxy = (win_size-1) / 2;
    x_min = x_cent - dxy; x_max = x_cent + dxy;
    y_min = y_cent - dxy; y_max = y_cent + dxy;
    x_mins = [x_mins, x_min]; x_maxs = [x_maxs, x_max];
    y_mins = [y_mins, y_min]; y_maxs = [y_maxs, y_max];
    windows(s, :) = [x_min, x_max, y_min, y_max];

    %% determine centroid in window
    ix_win = max(1, x_min):min(x_dim, x_max);
    iy_win = max(1, y_min):min(y_dim, y_max);
    window = E_img(iy_win, ix_win);
    spotcounts(s) = sum(window(:));
    window = window / spotcounts(s);
    x_loc = sum(sum(window, 1) .* linspace(ix_win(1), ix_win(end), ix_win(end)-ix_win(1)+1));
    y_loc = sum(sum(window, 2).' .* linspace(iy_win(1), iy_win(end), iy_win(end)-iy_win(1)+1));
    spots = [spots; [x_loc, y_loc]];
    win_weights(iy_win - y_cent + dxy + 1, ix_win - x_cent + dxy + 1, s) = window;
end

%% Prune spots that were detected too close to more prominent spots
[~, ixs] = sort(maxcorrvals, 'descend');
rm_ix = false(size(maxcorrvals));

for ix = ixs
    if ~rm_ix(ix)
        overlaps_bounds = abs(x_mins - x_mins(ix)) < (dxy+1) & ...  % If spot too close to another spot
            abs(y_mins - y_mins(ix)) < (dxy+1) & ...
            ~(abs(x_mins - x_mins(ix)) == 0 & (abs(y_mins - y_mins(ix)) == 0)); % Spot shouldn't remove itself
        overlaps_locs = abs(spots(:, 1) - spots(ix, 1)) < dxy & ...
            abs(spots(:, 2) - spots(ix, 2)) < dxy & ...
            ~(abs(spots(:, 1) - spots(ix, 1)) == 0) & ...
            ~(abs(spots(:, 2) - spots(ix, 2)) == 0);
        rm_ix(overlaps_bounds | overlaps_locs') = true;
    end
end

x_mins(rm_ix) = [];
y_mins(rm_ix) = [];
x_maxs(rm_ix) = [];
y_maxs(rm_ix) = [];
spots(rm_ix, :) = [];
windows(rm_ix, :) = [];
win_weights(:, :, rm_ix) = [];
maxcorrvals(rm_ix) = [];
spotcounts(rm_ix) = [];

%% plot detected spot windows
if plotDetectedSpots
    figure; imagesc(E_img); colorbar;
    set(gca, 'colorscale', 'log')
    hold on;
    for s=1:size(spots,1)
        x_min = x_mins(s); x_max = x_maxs(s);
        y_min = y_mins(s); y_max = y_maxs(s);
        x_loc = spots(s, 1); y_loc = spots(s, 2);
        % plot square for window
        plot(linspace(x_min, x_max, win_size), y_min*ones(win_size), 'r');
        plot(linspace(x_min, x_max, win_size), y_max*ones(win_size), 'r');
        plot(x_min*ones(win_size), linspace(y_min, y_max, win_size), 'r');
        plot(x_max*ones(win_size), linspace(y_min, y_max, win_size), 'r');
        % plot circle at centroid
        t=0:0.1:2*pi;
        x_circ = x_loc + 0.5*cos(t);
        y_circ = y_loc + 0.5*sin(t);
        plot(x_circ, y_circ, 'r' ,'LineWidth', 1);
    end
    hold off;
end
end
