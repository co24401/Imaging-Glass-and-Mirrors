filedir = '.\Data\big_mirror_collection\';
cscanfile =  '.\Data\big_mirror_collection\detector_scan_200x200.mat';
replicafile = '.\Data\replica\replica.mat';

paramfilename = 'params.mat';  % File to be saved as output

% Parameters for removing blur caused by galvo motion
dwell = 0.005; % seconds
lead_buffer = 0.0005; % seconds
trail_buffer = 0; % seconds
prf = 20E6; % Hz
first_count_only = false;

load(cscanfile, 'num_u', 'num_v', 'u_lims', 'v_lims')
num_pix = num_u*num_v;
num_bins = 8192;
num_pts = 100;

load(replicafile, 'mu_fit_LUT', 'pkbias_LUT')

savefigs = true;
savedata = true;

bin_radius = 30;
noise_gate_start = 1;
noise_gate_stop = 1000;
min_pk_cts_threshold = 0.7; %0.45; 

%min_pk_cts_thresholds = 0.45*ones(num_pts, 1);
%min_pk_cts_thresholds([1 2 4 20 37 39 40]) = 0.25;
%min_pk_cts_thresholds(21) = 0.18;

pkbins = zeros(num_pts, num_pix);
mubins = zeros(num_pts, num_pix);
muvars = zeros(num_pts, num_pix);
Evals = zeros(num_pts, num_pix);
Evars = zeros(num_pts, num_pix);
w2s = zeros(num_pts, num_pix);
S = false(num_pts, num_pix);

warning('off', 'MATLAB:nearlySingularMatrix')

%%

disp('Processing data collections.')

for ii = 1:num_pts 
    filename = ['spot_' num2str(ii)];
    
    disp(['Laser position ' num2str(ii)])
    
    [dataCube, ~, ~, ~] = processRawCounts( ...
        filename, filedir, savedata, savefigs, ...
        num_u, num_v, u_lims, v_lims, ...
        dwell, lead_buffer, trail_buffer, prf, num_bins, first_count_only);
    
    [pkbins(ii, :), mubins(ii, :), muvars(ii, :), Evals(ii, :), Evars(ii, :), ...
        w2s(ii, :), S(ii, :)] = datacubeStats(dataCube, ...
        noise_gate_start, noise_gate_stop, min_pk_cts_threshold, ...
        mu_fit_LUT, pkbias_LUT);

%     % Can alternatively define params like min_cts_threshold for each
%     frame independently.
%     [pkbins(ii, :), mubins(ii, :), muvars(ii, :), Evals(ii, :), Evars(ii, :), ...
%         w2s(ii, :), S(ii, :)] = datacubeStats(dataCube, ...
%         noise_gate_start, noise_gate_stop, min_pk_cts_thresholds(ii), ...
%         mu_fit_LUT, pkbias_LUT);
    
    close all
    
end

warning('on', 'MATLAB:nearlySingularMatrix')

save([filedir paramfilename], 'pkbins', ...
    'mubins', 'muvars', 'Evals', 'Evars', 'w2s', 'S', 'dwell', 'min_pk_cts_threshold')
