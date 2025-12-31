%RUN_EDCA_DEMO Configure and exercise an EDCA scenario with multiple users.
%   The script constructs a single-AP environment with several stations that
%   each generate traffic for the four EDCA access categories. It then plots
%   throughput and delay statistics to give a quick visual check of the MAC
%   behavior.

% Simulation configuration
cfg.numStations = 8;
cfg.totalSlots = 2e5; % keep reasonably small for quick runs
cfg.arrivalProb = [0.05 0.02 0.01 0.005]; % BE, BK, VI, VO arrival chance per slot

% Define standard-like EDCA parameters
% AIFSN and TXOP limits follow the IEEE 802.11 default EDCA parameter set
% (OFDM/high-throughput PHY defaults): AC_BK/AC_BE TXOP = 0, AC_VI = 3.008 ms,
% AC_VO = 1.504 ms. With 9 us slots, those map to 1,1,334,167 slots.
cfg.acParams = struct( ...
    'name',  {'AC_BE','AC_BK','AC_VI','AC_VO'}, ...
    'aifsn', {3, 7, 2, 2}, ...
    'cwMin', {15, 15, 7, 3}, ...
    'cwMax', {1023, 1023, 15, 7}, ...
    'txopSlots', {1, 1, 334, 167});

% Run simulation
results = edca_simulation(cfg);

% Display summary metrics
fprintf('Collisions: %d\n', results.collisions);
fprintf('Throughput per AC (Mbps)\n');
disp(array2table(results.throughputMbps', 'VariableNames', {cfg.acParams.name}));

fprintf('Average delay per AC (slots)\n');
disp(array2table(results.avgDelaySlots', 'VariableNames', {cfg.acParams.name}));

fprintf('Access success rate across all stations per AC (successes/attempts)\n');
disp(array2table(results.acSuccessRate', 'VariableNames', {cfg.acParams.name}));

fprintf('Access success rate per station (successes/attempts)\n');
disp(array2table(results.stationSuccessRate, 'VariableNames', {cfg.acParams.name}));

% Plot throughput per station for visualization
figure(1); clf;
bars = bar(results.stationThroughput, 'stacked');
title('Station Throughput by Access Category');
xlabel('Station index');
ylabel('Throughput (Mbps)');
legend(bars, {cfg.acParams.name}, 'Location', 'northoutside', 'Orientation', 'horizontal');
grid on;

% Plot medium usage over time (sampled subset)
figure(2); clf;
sampleRange = 1:min(500, numel(results.mediumUsage));
stem(sampleRange, results.mediumUsage(sampleRange), '.');
yticks([0 1 2]);
yticklabels({'Idle','Success','Collision'});
title('Medium Usage (0=Idle, 1=Success, 2=Collision)');
xlabel('Slot index');
ylabel('Channel state');
grid on;

% The figures are created for interactive inspection. In a headless
% environment the bar and stem data still help sanity-check the EDCA
% behavior (e.g., collisions vs. idle vs. busy proportions).
