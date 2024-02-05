% This is the actual Kalman filter
function [x,y,vx,vy,ax,ay,mm] = kfilter(posx,posy,post,order,Q,R,firstrun,missing)
% allocate memory for return values
n = length(posx);
x = zeros(n,1);
y = zeros(n,1);
vx = zeros(n,1);
vy = zeros(n,1);
ax = zeros(n,1);
ay = zeros(n,1);
mm = zeros(n,1);
% initialise return values and filtered state estimate
x(1) = posx(1);
y(1) = posy(1);
vx(1) = 0;
vy(1) = 0;
ax(1) = 0;
ay(1) = 0;
switch (order)
    case {0}
        cX = [posx(1) posy(1)]';
    case {1}
        cX = [posx(1) posy(1) 0 0]';
    case {2}
        cX = [posx(1) posy(1) 0 0 0 0]';
end
cP = 0.1*eye(2 + 2*order);
I = eye(2+2*order);
outlier = 0;
for k = 2:n
    % compute A and H from the time lag
    T = post(k) - post(k-1);
    switch(order)
        case {0}
            A = [1 0; ...
                0 1];
            H = [1 0; ...
                0 1];
        case {1}
            A = [1 0 T 0; ...
                0 1 0 T; ...
                0 0 1 0; ...
                0 0 0 1];
            H = [1 0 0 0; ...
                0 1 0 0]; ...
        case {2}
        A = [1 0 T 0 0.5*T*T 0; ...
            0 1 0 T 0 0.5*T*T; ...
            0 0 1 0 T 0; ...
            0 0 0 1 0 T; ...
            0 0 0 0 1 0; ...
            0 0 0 0 0 1;];
        H = [1 0 0 0 0 0; ...
            0 1 0 0 0 0];
    end
    if (firstrun & missing(k))
        % missing data, only predict
        % in the next EM steps they will be augmented with their MAP estimates
        pX = A * cX;
        pP = A * cP * A' + Q;
        % the equations are obtained by setting the Kalman gain to zero
        cX = pX;
        cP = pP;
    else % data not missing, predict and correct
        % prediction step
        pX = A * cX;
        pP = A * cP * A' + Q;
        % observation
        z = [posx(k); posy(k)];
        residual = z - H*pX;
        % validation gate for robustifying against extreme outliers
        chisq = [1000, 1000, 1000]; % just use a ridicilous threshold
        invS = (H*pP*H' + R)^-1; % measurement prediction covariance
        mahalanobis = residual' * invS * residual;
        mm(k) = mahalanobis;
        if ((outlier>5) | (mahalanobis < chisq(order+1)))
            % within validation gate -- perform correction step with
            % the new measurement (i.e. non-zero Kalman gain)
            K = pP * H' * invS;
            cX = pX + K*residual;
            cP = (I - K*H)*pP;
            outlier = 0;
        else
            % outlier -- ignore this measurement
            % the equations are obtained by setting the Kalman gain to zero
            cX = pX;
            cP = pP;
            outlier = outlier + 1;
        end
    end
    % save Kalman filtered states (cX)
    x(k) = cX(1);
    y(k) = cX(2);
    if (order > 0)
        vx(k) = cX(3);
        vy(k) = cX(4);
    end
    if (order > 1)
        ax(k) = cX(5);
        ay(k) = cX(6);
    end
end