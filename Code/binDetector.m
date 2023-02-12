% INPUTS:
%   pkbin_img   = Image of fitted peak bins - dim(y, x) (x, y = spatial coordinates)
%   pkvar_img   = Image of peak bin fit variances - dim(y, x)
%   windows     = sx4 matrix containing x_min, x_max, y_min, y_max locations for windows of s beams
%   t_0         = Time in receive window associated with zero time-of-flight
%   bin_width   = length of timing bin in seconds
%
% OUTPUTS: 
%   tofs        = tof for each spot (list of s entries)

function [tofs] = binDetector(pkbin_img, pkvar_img, windows, bin_width, t_0)
    x_dim = size(pkbin_img, 2);
    y_dim = size(pkbin_img, 1);

    % Used to correct biases in peak fitting
    LUT_dir = 'C:\Users\dolor\Remote-Flash\Experiments\Replica Collection 042821\replica.mat';
    load(LUT_dir, "mu_fit_LUT", "pkbias_LUT");
    
    %% loop through each spot
    tofs = [];
    for s=1:size(windows, 1) 
        %% extract data from window
        x_min = windows(s, 1); x_max = windows(s, 2); x_cent = (x_min+x_max)/2;
        y_min = windows(s, 3); y_max = windows(s, 4); y_cent = (y_min+y_max)/2;
        ix_win = max(1, x_min):min(x_dim, x_max);
        iy_win = max(1, y_min):min(y_dim, y_max);
        
        tof_vals = pkbin_img(iy_win, ix_win);
        tof_vars = pkvar_img(iy_win, ix_win);
        %% extract tof for each pixel in window
        % Using var^2 is strange.  Not sure if it is optimal.  But it applies
        % uncetainty weighting AND a beam-shape filter.  Using E for beam
        % shape filter may be more appropriate but it presents other issues
        % due to how the values are computed.  1/var is ~ equiv to E, but
        % possibly discourages noise-like returns even more.
        tof = bin_width*sum(tof_vals(:)./tof_vars(:).^2, 'omitnan')/sum(1./tof_vars(:).^2, 'omitnan') - t_0;
        tofs = [tofs, tof];
    end
end