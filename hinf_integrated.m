function [pos,vel,poshat, velhat, poshatinf, velhatinf, HinfGains, KalmanGains] = hinf_integrated(g, duration, dt, SteadyState)

coder.cinclude('hinf.h');

x = [0; 0]; % initial state vecto

% Initialize Kalman filter variables
xhat = x; % initial Kalman filter state estimate

% Initialize H-infinity filter variables
xhatinf = x; % initial H-infinity filter state estimate

% Initialize arrays for later plotting.
pos = [x(1)]; % true position array
vel = [x(2)]; % true velocity array
poshat = [xhat(1)]; % estimated position array (Kalman filter)
velhat = [xhat(2)]; % estimated velocity array (Kalman filter)
poshatinf = [xhatinf(1)]; % estimated position array (H-infinity)
velhatinf = [xhatinf(2)]; % estimated velocity array (H-infinity)
HinfGains = [0; 0]; % H-infinity filter gains
KalmanGains = [0; 0]; % Kalman filter gains

[pos,vel,poshat, velhat, poshatinf, velhatinf, HinfGains, KalmanGains] = coder.ceval('hinf', g, duration, dt, SteadyState);

end

