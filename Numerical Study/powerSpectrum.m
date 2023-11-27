T = readtable('turbulence1.csv');	                                        % Load in the saved data
t = T{:,1};
f = T{:,2};

[pow, f] = pspectrum(f);                                                  % Create the powerspectrum

set(groot,'defaultAxesTickLabelInterpreter','latex')

figure(1), clf(1)
loglog(f/pi, pow)                                                         % Create a loglog plot
xlim([0.001 1])                                                           % Set reasonable x-limits
xlabel('Frequency','Interpreter','latex'), ylabel('Spectral power','Interpreter','latex')
