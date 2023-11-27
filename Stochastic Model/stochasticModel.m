%% Stochastic model for large scale convection roll
clear, clc, clf, format shortG

L = 0.2476;                                                                 % m
nu = 6.69e-3/1e4;                                                           % cm^2/s --> m^2/s
Re = 3700;                                                                  % -
D_d = 3.5e-5;                                                               % K^2/s, subscript d stands for delta
D_td = 2.5e-5;                                                              % rad^2/s^3, subscript td stands for theta_dot

tau_d = L^2/(18*nu*sqrt(Re));                                               % tau_delta
tau_td = L^2/(2*nu*Re);                                                     % tau_theta_dot
delta0 = 0.25;                                                              % K, constant value
delta_0 = delta0;                                                           % Changes in the iteration
thetad_0 = 0;                                                               % rad/s
theta_0 = 0;                                                                % rad

dt = 1;                                                                     % s, timestep
t_end = 60*60*24;                                                           % s, 1 day simulation time
t_plot = dt:dt:t_end;                                                       % Create a list for plotting

f_d = normrnd(0,sqrt(D_d/dt));                                              % Initial value of the delta stochastic term
f_td = normrnd(0,sqrt(D_td/dt));                                            % Initial value of the theta_dot stochastic term

delta = zeros(1,length(t_plot));                                            % Initialization
thetad = zeros(1,length(t_plot));                                           % Initialization
theta = zeros(1,length(t_plot));                                            % Initialization

for i = 1:length(t_plot)
    delta(i) = dt*((delta_0/tau_d)*(1-sqrt(delta_0/delta0))+f_d)+delta_0;   % Time integration of delta
    
    delta_0 = delta(i);                                                     % Updating old value
    f_d = normrnd(0,sqrt(D_d/dt));                                          % Calculating new stochastic term         
end

figure(1)                                                                   % Plotting
plot(t_plot,delta)
xlabel('Time [s]'), ylabel('\delta [K]')

for j = 1:length(t_plot)
    theta(j) = theta_0 + thetad_0*dt;                                       % Time integration of theta
    thetad(j) = thetad_0 + (f_td-(thetad_0*delta(j)/(tau_td*delta0)))*dt;   % Time integration of theta_dot

    theta_0 = theta(j);                                                     % Updating old value
    thetad_0 = thetad(j);                                                   % Updating old value
    f_td = normrnd(0,sqrt(D_td/dt));                                        % Calculating new stochastic term 
end

figure(2)                                                                   % Plotting
plot(t_plot,theta)
xlabel('Time [s]'), ylabel('\theta_0/2\pi')

figure(3), clf(3), hold on                                                  % Plotting
T_0 = 293;                                                                  % K, assumption
T = T_0 + delta.*cos(0 - theta);                                            % Taking theta = 0
plot(t_plot,T)
xlabel('Time [s]'), ylabel('T [K]')
hold off
