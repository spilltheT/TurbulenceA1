%% Lagrangian description of flow
clear, clc, clf, format longE

%% Exercise b)

% k = 257;                                                                    % Amount of points on the circle
% th = linspace(0,2*pi,k);                                                        
% x = cos(th);                                                                % Creating the x and y-coordinates
% y = sin(th);
% 
% x_new = zeros(1,length(x));                                                 % Initializing vector for new position x-values
% y_new = zeros(1,length(y));                                                 % Initializing vector for new position y-values
% dt = 0.001;                                                                  % Time step
% t = dt;                                                                     % Variable used in the iteration
% t_end = [0 0.1 1 10];                                                       % End time
% 
% u_x = @(y,t) 2*y*cos(t);                                                    % Velocity in x-direction
% u_y = @(x,t) 4*x*sin(2*t);                                                  % Velocity in y-direction
% 
% for j = 1:length(t_end)                                                       % Iterating over the different end times
%     for k = 0:dt:t_end(j)
%         for l = 1:length(x)
%             x_new(l) = x(l) + integral(@(t) u_x(y(l),t),t-dt,t);            % Determining new position by integrating velocity
%             y_new(l) = y(l) + integral(@(t) u_y(x(l),t),t-dt,t);
% 
%             x(l) = x_new(l);                                                % Replace old values for next iteration
%             y(l) = y_new(l);
%         end
%         t = t + dt;                                                         % Add timestep to current time for next iteration
%     end
%     figure(j), hold on                                                      % Plot the newest positions of the dye (at t_end)
%     set(gcf, 'Position',  [100, 100, 600, 600])
%     plot(x_new,y_new), daspect([1 1 1]),  grid on
%     xlabel('{\it x}'), ylabel('{\it y}'), title(['Dye in a stream function at{\it t} = ' num2str(t_end(j))])
%     hold off
% end

%% Exercise c)

k = 257;
th = linspace(0,2*pi,k);
x = cos(th);
y = sin(th);

x_0 = x;
y_0 = y;

x_new = zeros(length(x),1);
y_new = zeros(length(y),1);
dt = 0.001;
t = dt;
t_end = pi;

u = @(y,t) 2*y*cos(t);
v = @(x,t) 4*x*sin(2*t);

for k = 0:dt:t_end
    for l = 1:length(x)
        x_new(l) = x(l) + integral(@(t) u(y(l),t),t-dt,t);
        y_new(l) = y(l) + integral(@(t) v(x(l),t),t-dt,t);

        x(l) = x_new(l);
        y(l) = y_new(l);
    end
    t = t + dt;
end                                                                         % All of the above the same as b) to initialize the flow field
                                                                            % Only difference is one less for-loop as only one t_end is required

tau = [0.1 1 10 100];                                                       % Tau vector

for j = 1:length(tau)

    t = dt;                                                                 % Variable used in the iteration
                                                                            % All below is initialized in the for-loop to ensure they are reset to their original values once tau changes
    x_p = zeros(length(x),1);                                               % Subscript p for particle (of the dye)
    y_p = zeros(length(x),1);
    u_p = zeros(length(x),1);                                               % Current particle velocity
    v_p = zeros(length(x),1);
    u_po = zeros(length(x),1);                                              % Old particle velocity
    v_po = zeros(length(x),1);

    for i = 1:length(x)
        u_po(i) = u(y_0(i),0);                                              % Initial particle velocities (being equal to the flow velocity)
        v_po(i) = v(x_0(i),0);
    end

    k = 257;
    th = linspace(0,2*pi,k);
    x_po = cos(th);                                                         % Initial particle position on the unit circle 
    y_po = sin(th);

    for k = 0:dt:t_end
        for l = 1:length(x)
            x_p(l) = x_po(l) + dt*u_po(l);                                  % Standard time integration, using the old velocity value       
            y_p(l) = y_po(l) + dt*v_po(l);
            u_p(l) = u_po(l) + (dt/tau(j))*(u(y_p(l),t)-u_po(l));
            v_p(l) = v_po(l) + (dt/tau(j))*(v(x_p(l),t)-v_po(l));

            x_po(l) = x_p(l);                                               % Updating old values
            y_po(l) = y_p(l);
            u_po(l) = u_p(l);
            v_po(l) = v_p(l);
        end
        t = t + dt;
    end
    figure(1), hold on                                                      % Plotting
    plot(x_p,y_p), grid on%, axis([-20 20 -10 10])
    xlabel('{\it x}'), ylabel('{\it y}'), title(['Dye in a stream function at{\it t} = ' num2str(t_end)])
end

plot(x_new,y_new,'--'), legend('{\it\tau} = 0.1', '{\it\tau} = 1', '{\it\tau} = 10', '{\it\tau} = 100', 'Original flow')
hold off
