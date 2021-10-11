
%% time window
T_start=-20;
T_end=20;
dt=0.01;

%% target distribution
    theta_target=[0.6,7,0.4^2,5,0.3^2];
    density_target=mixture_density(theta_target,T_start:dt:T_end);
    % ideal case, know the density
    
%% initial distribution, parameterized by mixture family (theta)
theta_initial=[0.3,-3,0.5^2,-5,0.4^2];

%% W2obj_minimization

% 5 methods:
% WGD Wasserstein GD
% FRGD Fisher-Rao GD
% GD original GD
% pGD preconditioned GD, with G(1,1)=30;
% WGDm modified WGD

method='WGD';
mixture_W2obj_minimization( T_start, T_end, dt, theta_initial, density_target,method )
% will output the history of the iterations
% and output the initial, true, and the converged solution

