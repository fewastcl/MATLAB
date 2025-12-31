%RUN_EDCA_SCALING Sweep user count and plot access success rates.
%   This script runs the EDCA simulator for user counts 1..50, records the
%   per-access-category access success rate (successes/attempts averaged
%   across stations), and plots how the probability of a successful access
%   changes as more users contend for the channel.

userCounts = 1:50;
numSweeps = numel(userCounts);
acSuccessRate = zeros(numSweeps, 4); % rows: user count, cols: BE BK VI VO

baseConfig.totalSlots = 1e5;
baseConfig.arrivalProb = [0.05 0.02 0.01 0.005]; % BE, BK, VI, VO
baseConfig.acParams = struct( ...
    'name',  {'AC_BE','AC_BK','AC_VI','AC_VO'}, ...
    'aifsn', {3, 7, 2, 2}, ...
    'cwMin', {15, 15, 7, 3}, ...
    'cwMax', {1023, 1023, 15, 7}, ...
    'txopSlots', {8, 4, 4, 2});

for idx = 1:numSweeps
    cfg = baseConfig;
    cfg.numStations = userCounts(idx);
    results = edca_simulation(cfg);

    % Average per-station success rate for each AC (omit NaNs for stations
    % that never attempted to transmit in an AC).
    acSuccessRate(idx, :) = mean(results.stationSuccessRate, 1, 'omitnan');
end

figure(3); clf;
plot(userCounts, acSuccessRate(:, 1), '-o', ...
     userCounts, acSuccessRate(:, 2), '-s', ...
     userCounts, acSuccessRate(:, 3), '-^', ...
     userCounts, acSuccessRate(:, 4), '-d');
xlabel('Number of stations');
ylabel('Access success rate (successes / attempts)');
title('EDCA Access Success Rate vs. User Count');
legend({'AC\\_BE','AC\\_BK','AC\\_VI','AC\\_VO'}, 'Location', 'southwest');
grid on;
