function [X_0,X_0_record] = mirrorSourceFit(XoB, toB, X_0_init, t, epsilon, max_iter, plot_trajectory)
% mirrorSourceFit.m
% Connor Henley
% 3/4/2022
%
% Robustly localizaes the position of the mirrored source using a 
% RANSAC algorithm.
%
% Inputs:
% XoB = Set of out of beam (two-bounce) points (3xN)
% toB = Times of flight associated with out of beam points (1xN)
% X_0_init = Initial guess for mirrored source position (3x1)
% t = Damping constant for Newton step (1 = undamped).
% epsilon = Stopping criterion.
% max_iter = Maximum number of Newton steps.
% plot_trajectory = If true, plot trajectory of mirror source estimate vs.
%   number of iterations.

% Outputs:
% X_0 = Estimated mirrored source positin (3x1)
% X_0_record = History of mirrored source position estimates (3x[max_iter+1])

c = 299792458; % speed of light (m/s)

X_0 = X_0_init;
r_xc_oB = vecnorm(XoB, 2, 1);

X_0_record = zeros(3, max_iter+1);
X_0_record(:, 1) = X_0;

lambda_record = zeros(size(1, max_iter));

for ii = 1:max_iter
    
    dx = X_0(1) - XoB(1,:);
    dy = X_0(2) - XoB(2,:);
    dz = X_0(3) - XoB(3,:);
    dt = c*toB - r_xc_oB;
    
    g = dx.^2 + dy.^2 + dz.^2 - dt.^2;
    
    grad_array = 2*g.*[dx; dy; dz]; % Gradient
    
    % Hessian terms
    h11 = reshape(dx.^2, 1, 1, []);
    h12 = reshape(dx.*dy, 1, 1, []);
    h13 = reshape(dx.*dz, 1, 1, []);
    h22 = reshape(dy.^2, 1, 1, []);
    h23 = reshape(dy.*dz, 1, 1, []);
    h33 = reshape(dz.^2, 1, 1, []);
    
    % Hessian
    H_array = 4*[h11 h12 h13; h12 h22 h23; h13 h23 h33] + ...
        2*reshape(g, 1, 1, []).*[1 0 0 ; 0 1 0 ; 0 0 1];
    
    grad_f = sum(grad_array, 2);
    H_f = sum(H_array, 3);
    
    delta = -inv(H_f)*grad_f;
    lambda = -0.5*grad_f'*delta;
    lambda_record(ii) = lambda;
    
    X_0 = X_0 + t*delta;
    X_0_record(:, ii+1) = X_0;
    
    e = eig(H_f)';
    if sum(e<0) > 0
        disp(['Warning, H was not positive definite at step' num2str(ii)])
    end
    
    if abs(lambda) < epsilon
        break
    end
    
end

disp(['Finished after ' num2str(ii) ' iterations.'])

if plot_trajectory

    figure; scatter3(X_0_record(1,1:(ii+1)), X_0_record(2,1:(ii+1)), X_0_record(3,1:(ii+1)), 10, 1:(ii+1), 'filled')
    hold on; scatter3(X_0_record(1,1), X_0_record(2,1), X_0_record(3,1), 80, 'kp' ,'filled')
    hold on; scatter3(X_0_record(1,ii+1), X_0_record(2,ii+1), X_0_record(3,ii+1), 80, 'rp' ,'filled')
    legend({'Guesses', 'Init.', 'Final'})
    colorbar;
    
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    
    axis equal
    grid on
    
    view(-165, -60)
end
end

