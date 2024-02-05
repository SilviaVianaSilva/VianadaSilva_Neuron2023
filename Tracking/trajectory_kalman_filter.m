function [t,x,y,vx,vy,ax,ay] = trajectory_kalman_filter(posx,posy,post,order,Q,R)
%
% Kalman filter for obtaining an appriximate Bayesian MAP estimate  to the rat's
% trajectory from raw video tracker data.
%
% The trajectory model is the set of Taylor series expansions:
%
%    x(t) = x(a) + x'(a)*(t-a) + 1/2 * x''(a)*(t-a)^2 + o(h)
%    x'(t) = x'(a) + x''(a)*(t-a) + o(h)
%    x''(t) = x''(a) + o(h)
%
% The prediction step vill approximate the state variables [x(t), x'(t), x''(t)]
% by expanding around a = t-1. The correction step will use the observation
% at time t and the estimated state to create a corrected ("filtered") estimate
% of the state at time t. The corresponding state space model is:
%
%    state:          x(t+1) = A * x(t) + w, with w ~ N(0,Q)
%    observation:    z(t) = H * x(t) + v,   with v ~ N(0,R)
%
%
% Mandatory input arguments to the routine:
%
%    posx, posy, post -- tracker raw data and time stamps
%
%    order -- order of approximation (Taylor series) to the trajectory
%       0: ignore all derivatives (locally constant approximation)
%       1: include first derivatile (locally linear approximation)
%       2: include first and 2nd derivative (locally quadratic approximation)
%
%       NB! It looks like the quadratic approximation is the better method, at
%       least for recordings where the rat swims in circle(s). Note that circular
%       velocity and acceleration are not included in the state model, we might
%       get a better filter for these data by including these parameters.
%
% Optional arguments:
%
%    Q -- a priori specified system noise covariance matrix, rank (2+2*order) x (2+2*order)
%       default value: 0.1 * eye(2 + 2*order)
%
%    R -- a priori specified measurement noise covariance matrix, rank 2 x 2
%       default value: eye(2)
%
% Output arguments from the Kalman filter:
%
%    t(k) -- time stamps for filter output
%
%    [x(k), y(k)] -- approximately the most probable location at time t(k) conditional
%       on the observed [posx,posy] in the period [t(1),t(k)], the a priori specified
%       Q and R matrices, and the prior state variable.
%
%    [vx(k), vy(k)] -- the most probable velocity at time t(k) conditional on the
%       observed [posx,posy] in the period [t(1),t(k)] and the a priori specified Q and
%       R matrices, and the prior state variable.
%
%    [ax(k), ay(k)] -- the most probable acceleration at time t(k) conditional on the
%       observed [posx,posy] in the period [t(1),t(k)], the a priori specified Q and
%       R matrices, and the prior state variable.
%
%    Missing positions indicated by the value NaN for posx and posy are fitted using
%    the EM algorithm of Dempster et al. (1977).
%
%
% Copyright (C) 2004 Sturla Molden
% Centre for the Biology of Memory
% Norwegian University of Science and Technology
%
%

% Default values. These are manually tuned for good performance
if (nargin < 5)
    Q = 0.1 * eye(2 + 2*order);
end
if (nargin < 6)
    R = eye(2);
end
% Chop off any missing data in the start of session
n = length(posx);
lastmissing = min(find(isfinite(posx))) - 1;
posx = posx(lastmissing+1:end);
posy = posy(lastmissing+1:end);
post = post(lastmissing+1:end);
% Run Kalman filter on the remaining samples
missing = zeros(n,1);
missing(isnan(posx)) = 1;
if (sum(missing))
    % Missing data, use EM algorithm to get MAP estimators (Dempster et al. 1977)
    % Find the missing data we want to augment, initially augment
    % with Kalman predicted positions, then augment with Kalman filtered
    % positions. Iterate ten times to allow the augmented data to converge.
    missing_index = find(missing);
    [x,y,vx,vy,ax,ay] = kfilter(posx,posy,post,order,Q,R,1,missing);
    for i = 1:10
        posx(missing_index) = x(missing_index);
        posy(missing_index) = y(missing_index);
        [x,y,vx,vy,ax,ay] = kfilter(posx,posy,post,order,Q,R,0);
    end
else
    % No missing data, get the MAP estimates from a single pass
    [x,y,vx,vy,ax,ay,mm] = kfilter(posx,posy,post,order,Q,R,0);
end
t = post;