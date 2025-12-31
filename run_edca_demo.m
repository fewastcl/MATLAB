%RUN_EDCA_DEMO Configure and exercise an EDCA scenario with multiple users.
%   The script constructs a single-AP environment with several stations that
%   each generate traffic for the four EDCA access categories. It then plots
%   throughput and delay statistics to give a quick visual check of the MAC
%   behavior.

% Simulation configuration
cfg.numStations = 8;
cfg.totalSlots = 2e5; % keep reasonably small for quick runs
cfg.arrivalProb = [0.02 0.05 0.01 0.005]; % BK, BE, VI, VO arrival chance per slot
cfg.slotTime = 9e-6; % slot duration (s)
cfg.phyRateMbps = [600 600 1200 1200]; % per-AC PHY data rates (example MCS/NSS)
cfg.phyPreambleUs = 100; % HE preamble duration (us)
cfg.macHeaderBits = 272; % MAC overhead per PPDU
cfg.collisionDurationUs = 100; % wasted medium time on collision

% Define standard-like EDCA parameters
% AIFSN and TXOP limits follow the IEEE 802.11 HE default EDCA parameter set:
%   AC_BK/AC_BE TXOP = 2528 us, AC_VI = 4096 us, AC_VO = 2080 us.
% TXOP values are expressed directly in microseconds to match the standard.
cfg.acParams = struct( ...
    'name',  {'AC_BK','AC_BE','AC_VI','AC_VO'}, ...
    'aifsn', {7, 3, 2, 2}, ...
    'cwMin', {15, 15, 7, 3}, ...
    'cwMax', {1023, 1023, 15, 7}, ...
    'txopUs', {2528, 2528, 4096, 2080});

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
legendLabels = {'AC_{BK}','AC_{BE}','AC_{VI}','AC_{VO}'};
legend(bars, legendLabels, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex');
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
