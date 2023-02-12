function [pD, pS1, pS2, diffuse_first] = computeGeometryMirrorDisambiguation(detections, s, refl_cal)
% computeGeometryMirrorDisambiguation.m
% Connor Henley
% 3/11/2022
%
% Compute 3D points from detections.

% Inputs:
% detections = Structure that information about spots detected in
% the frame, such as time and angle of arrival for each spot.
% s = Baseline seperation between source (L) and camera aperture (origin).
% refl_cal = Unitless constant that enables computation of reflectance from
% distance to and intensity received from a spot.
%
% Outputs:
% pD = Structure that holds information about point on diffusely reflecting
% surface (D)
% pS1 = Structure that holds information about point on specular surface
% that is illuminated with the laser (S1)
% pS2 = Structure that holds information about points on specular surfaces
% that are observed by the detector (S2).  Also
% diffuse_first = Flag that is true if diffusely refelcting surface was
% illuminated first, and false otherwise.

c = 299792458; % speed of light (m/s)
L = [s; 0; 0]; % source position
in_beam_tolerance = 0.0001; %0.01; %0.001; % Maximum cosine distance of point from beam ray
%   for point to classified as lying along laser beam path

% Load detection params to de-clutter code
tofs = detections.tofs;
thetas = detections.thetas;
phis = detections.phis;
spotcounts = detections.spotcounts;
ltheta = detections.ltheta;
lphi = detections.lphi;

pD = struct('pos', [], 'r', [], 'theta', [], 'phi', [], ...
    'tof', [], 'ix', [], 'rl', [], 'counts', [], 'reflectance', [], ...
    'behind_window', []);
pS1 = struct('pos', [], 'n', [], 'rl', []);
pS2 = struct('pos', [], 'n', [], 'r', [], 'thetas', [], ...
    'phis', [], 'tofs', [], 'counts', [], 'reflectance', []);
diffuse_first = 0;

if ~isempty(tofs)
    
    % Naively compute positions of all spots assuming that they've followed
    % one-bounce trajectories.
    r_xc = 0.5*(c^2 * tofs.^2 - s^2)./( c*tofs - s*cos(thetas) );
    X = r_xc.*[cos(thetas); sin(thetas).*sin(phis); sin(thetas).*cos(phis)];
    dir_point_L = [cos(ltheta); sin(ltheta)*sin(lphi); sin(ltheta)*cos(lphi)]; % Beam vector
    
    dir_XL_est = (X - L) ./ vecnorm(X - L, 2, 1); % Directions from laser position to naively computed spot positions.
    sepFromBeam = (1 - dir_point_L'*dir_XL_est); % Cosine distances from (apparent angle to naive pos's from L) and beam vector
    inBeam = sepFromBeam < in_beam_tolerance; % Spots within cos distance tol. of beam vector
    
    if sum(inBeam) > 1 % Window disambiguation step triggered if multiple detections along beam
        diffuse_first = 0;
        thru_window = 1;
        if sum(inBeam)>2
            [~, bix] = mink(sepFromBeam, 2); % If g.t. two spots near beam, choose two nearest
        else
            bix = find(inBeam);
        end
        
        if sum(inBeam) < length(inBeam)
            [dtof, dix] = min(tofs./(~inBeam), [], 'omitnan'); % Index of "true"
        % spot that is inside room.  This method of choosing dix will fail
        % if D is inBeam but further away from beam than two chosen bix
        % points.
        else
            dix = []; 
            dtof = Inf; % This will only ever 
            % be triggered if all detected spots are in beam.
        end
        
        if tofs(bix(1)) > dtof && tofs(bix(2)) < dtof
            s2ix = bix(1); % Index of spot observed on window that is also illuminated directly
            doutix = bix(2); % Index of spot behind window
        elseif tofs(bix(1)) < dtof && tofs(bix(2)) > dtof
            s2ix = bix(2);
            doutix = bix(1);
        elseif tofs(bix(1)) > dtof && tofs(bix(2)) > dtof
            refs = spotcounts(bix).*r_xc(bix).^2/refl_cal;
            if refs(2) > refs(1)
                s2ix = bix(1); % 3B reflection should be dimmer than 1B
                doutix = bix(2); % reflection off of point behind window
            else
                s2ix = bix(2); % 3B reflection should be dimmer than 1B
                doutix = bix(1); % reflection off of point behind window
            end
        else
            % If neither spot arrives after light from D,
            % then neither can be S2.  So we treat each as behind-window
            % detections.  This might get triggered if we illuminate an
            % edge.
            diffuse_first = 1;
            if tofs(bix(1)) < tofs(bix(2))
                dix = bix(1);
                dtof = tofs(bix(1));
                doutix = bix(2); 
            else
                dix = bix(2);
                dtof = tofs(bix(2));
                doutix = bix(1);
            end
        end
        
        if isempty(dix)
            diffuse_first = 1;
            dix = doutix(1);
            six = false(size(tofs));
        else
            six = true(size(tofs));
            six([dix doutix]) = false;
        end
    elseif sum(inBeam) == 1  % Only one detection along beam, use regular spot classification
        thru_window = 0;
        bix = find(inBeam);
        %btof = tofs(bix);
        [dtof, dix] = min(tofs, [], 'omitnan');
        six = tofs > dtof;
        if bix == dix
            diffuse_first = 1;
        else % Note, have already determined that dtof < btof at this point, dtof is min and dtof ~= btof
            diffuse_first = 0;
            bref = spotcounts(bix).*r_xc(bix).^2/refl_cal;
            dref = spotcounts(dix).*r_xc(bix).^2/refl_cal;
            if bref < dref
                s2ix = find(inBeam, 1);
            else
                thru_window = 1;
                dix = bix;
                six(bix) = false;
                diffuse_first = 1;
            end
        end
    else % No detections in beam means we have no way to compute depths
        return
    end
    
    dtheta = thetas(dix);
    dphi = phis(dix);
    
    % I hacked this in to deal with the case where specular surface was
    % illuminated first but we only see the diffuse spot.  There's almost
    % certainly a better way to deal with this that also considers case
    % where there are no spots in beam, but potentially many spots.
    if ~diffuse_first && length(tofs) == 1 % ~any(six)%
        % Store information about D
        pD.theta = dtheta;
        pD.phi = dphi;
        pD.tof = dtof;
        pD.ix = dix;
        pD.counts = spotcounts(dix);
        pD.behind_window = false;
        
        return;
    end
    
    stofs = tofs(six); % S spot ToFs
    sthetas = thetas(six); % S spot azimuths
    sphis = phis(six); % S spot elevations
    sspotcounts = spotcounts(six); % S spot counts
    
    cos_deltas = cos(dtheta)*cos(sthetas) + ...
        sin(dtheta)*sin(sthetas).*cos(dphi-sphis);
    
    dt_sd = stofs - dtof;
    
    if diffuse_first
        r_dc = r_xc(dix);
        D = X(:, dix);
    else
        r_dMc_s2 = r_xc(s2ix);
        r_dc = r_dMc_s2 - c*(tofs(s2ix) - dtof);
        D = r_dc*[cos(dtheta); ...
            sin(dtheta)*sin(dphi); ...
            sin(dtheta)*cos(dphi)];
    end
    
    refD = spotcounts(dix)*r_dc^2/refl_cal; % reflectance
    r_dl = norm(D - L);
    
    % Compute points on specular surface
    r_sc = (c/2)*( dt_sd.*(dt_sd + 2*r_dc/c) ./ ...
        ( dt_sd + (1-cos_deltas)*r_dc/c ));
    
    S = [r_sc.*cos(sthetas); ...
        r_sc.*sin(sphis).*sin(sthetas); ...
        r_sc.*cos(sphis).*sin(sthetas)];
    
    r_dMc = c*dt_sd + r_dc;
    refS = sspotcounts.*r_dMc.^2/refl_cal;
    
    % Compute surface normals
    if isempty(S)
        nS = [];
    else
        nS = -S ./ r_sc - (S - D) ./ (c*dt_sd + r_dc - r_sc);
        nS = nS ./ vecnorm(nS, 2, 1);
    end
    
    % Store information about D
    pD.pos = D;
    pD.r = r_dc;
    pD.theta = dtheta;
    pD.phi = dphi;
    pD.tof = dtof;
    pD.ix = dix;
    pD.rl = r_dl;
    pD.counts = spotcounts(dix);
    pD.reflectance = refD;
    if thru_window && diffuse_first
        pD.behind_window = true;
    else
        pD.behind_window = false;
    end
    
    % Store information about S points observed by detector
    pS2.pos = S;
    pS2.r = r_sc;
    pS2.n = nS;
    pS2.thetas = sthetas;
    pS2.phis = sphis;
    pS2.tofs = stofs;
    pS2.counts = sspotcounts;
    pS2.reflectance = refS;  
    
    if ~diffuse_first
        % Compute point on specular surface illuminated by laser
        r_dMl = c*dtof - r_dc; % Range from laser to mirror image spot
        
        dir_DL_est = (D - L) / r_dl; % Direction from laser to estimated D position
        cos_deltaL = dir_DL_est'*dir_point_L;
        
        r_s1l = 0.5*(r_dMl.^2 - r_dl.^2)./(r_dMl - r_dl*cos_deltaL);
        S1 = L+ r_s1l*[cos(ltheta); ...
            sin(ltheta)*sin(lphi); ...
            sin(ltheta)*cos(lphi)];
        
        % Compute surface normal
        nS1 = (L-S1) / r_s1l + (D - S1) / norm(D-S1);
        nS1 = nS1 / norm(nS1);
        
        pS1.pos = S1;
        pS1.n = nS1;
        pS1.rl = r_s1l;
    end
    
    if thru_window && ~diffuse_first
        douttheta = thetas(doutix);
        doutphi = phis(doutix);
        r_doutc = r_xc(doutix);
        Dout = r_doutc.*[cos(douttheta); ...
            sin(douttheta)*sin(doutphi); ...
            sin(douttheta)*cos(doutphi)];
        refDout = spotcounts(doutix).*r_doutc.^2/refl_cal; % reflectance
        r_doutl = vecnorm(Dout - L, 2, 1);
        
        pD.pos = [pD.pos Dout];
        pD.r = [pD.r r_doutc];
        pD.theta = [pD.theta douttheta];
        pD.phi = [pD.phi doutphi];
        pD.tof = [pD.tof tofs(doutix)];
        pD.ix = [pD.ix doutix];
        pD.rl = [pD.rl r_doutl];
        pD.counts = [pD.counts spotcounts(doutix)];
        pD.reflectance = [pD.reflectance refDout];
        pD.behind_window = [pD.behind_window true(size(doutix))];
    end
end

end

