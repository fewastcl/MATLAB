%RUN_EDCA_SCALING Sweep user count and plot access success rates.
%   This script runs the EDCA simulator for user counts 1..50, records the
%   per-access-category access success rate (successes/attempts averaged
%   across stations), and plots how the probability of a successful access
%   changes as more users contend for the channel.

userCounts = 1:50;
numSweeps = numel(userCounts);
acSuccessRate = zeros(numSweeps, 4); % rows: user count, cols: BE BK VI VO
totalThroughput = zeros(numSweeps, 1);
acThroughput = zeros(numSweeps, 4); % Mbps per AC

baseConfig.totalSlots = 1e5;
baseConfig.arrivalProb = [0.05 0.02 0.01 0.005]; % BE, BK, VI, VO
baseConfig.slotTime = 32e-6; % HE/EHT TXOP limits are multiples of 32 us
baseConfig.acParams = struct( ...
    'name',  {'AC_BE','AC_BK','AC_VI','AC_VO'}, ...
    'aifsn', {3, 7, 2, 2}, ...
    'cwMin', {15, 15, 7, 3}, ...
    'cwMax', {1023, 1023, 15, 7}, ...
    'txopSlots', {79, 79, 128, 65}); % HE TXOP limits (us) mapped to 32 us slots

for idx = 1:numSweeps
    cfg = baseConfig;
    cfg.numStations = userCounts(idx);
    results = edca_simulation(cfg);

    % Overall (all stations) access success rate per AC.
    acSuccessRate(idx, :) = results.acSuccessRate;

    % Throughput across all ACs for all users.
    acThroughput(idx, :) = results.throughputMbps';
    totalThroughput(idx) = sum(results.throughputMbps);
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

figure(4); clf;
plot(userCounts, totalThroughput, '-o', ...
     userCounts, acThroughput(:, 1), '-s', ...
     userCounts, acThroughput(:, 2), '-^', ...
     userCounts, acThroughput(:, 3), '-d', ...
     userCounts, acThroughput(:, 4), '-x');
xlabel('Number of stations');
ylabel('Aggregate throughput (Mbps)');
title('Throughput vs. User Count (Total and Per AC)');
legend({'Total','AC\\_BE','AC\\_BK','AC\\_VI','AC\\_VO'}, 'Location', 'northwest');
grid on;
