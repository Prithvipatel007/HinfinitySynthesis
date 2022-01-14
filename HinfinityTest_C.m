% Plot the results
close all; % Close all open figures

g = 0.01;
duration = 20;
dt = 0.1;
SteadyState = 0;


[pos,vel,poshat, velhat, poshatinf, velhatinf, HinfGains, KalmanGains] = Hinfinity_mex(g, duration, dt, SteadyState);
fprintf('Interpretation from Matlab Coder auto-gen C code');

t = 0 : dt : duration; % Create a time array

% Plot the position estimation error
% (Kalman filter = red line, H-infinity filter = green line)
figure;
plot(t,pos-poshat,'r', t,pos-poshatinf,'b--');
set(gca,'FontSize',12); set(gcf,'Color','White');
grid;
xlabel('Time (sec)');
ylabel('Position Error (feet)');
title('Position Estimation Error');
legend('Kalman filter', 'H_{\infty} filter');

% Plot the velocity estimation error
% (Kalman filter = red line, H-infinity filter = green line)
figure;
plot(t,vel-velhat,'r', t,vel-velhatinf,'b--');
set(gca,'FontSize',12); set(gcf,'Color','White');
grid;
xlabel('Time (sec)');
ylabel('Velocity Error (feet)');
title('Velocity Estimation Error');
legend('Kalman filter', 'H_{\infty} filter');

% Plot the Kalman filter gain matrix
figure;
plot(t,KalmanGains(1,:),'r', t,KalmanGains(2,:),'b--');
set(gca,'FontSize',12); set(gcf,'Color','White');
grid;
xlabel('Time (sec)');
title('Kalman Gains');

% Plot the H-infinity filter gain matrix
figure;
plot(t,HinfGains(1,:),'r', t,HinfGains(2,:),'b--');
set(gca,'FontSize',12); set(gcf,'Color','White');
grid;
xlabel('Time (sec)');
title('H-Infinity Gains');