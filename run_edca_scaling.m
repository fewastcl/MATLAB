%RUN_EDCA_SCALING Sweep user count and plot access success rates.
%   This script runs the EDCA simulator for user counts 1..50, records the
%   per-access-category access success rate (successes/attempts averaged
%   across stations), and plots how the probability of a successful access
%   changes as more users contend for the channel.

userCounts = 1:50;
numSweeps = numel(userCounts);
acSuccessRate = zeros(numSweeps, 4); % rows: user count, cols: BK BE VI VO
totalThroughput = zeros(numSweeps, 1);
acThroughput = zeros(numSweeps, 4); % Mbps per AC
busySuccess = zeros(numSweeps, 1);
busyCollision = zeros(numSweeps, 1);

% PHY configuration (set per-AC MCS/NSS/GI; rates are derived inside edca_simulation)
mcsPerAc = [3 3 7 9];        % UHR-MCS indices (0-based) for BK, BE, VI, VO
nssPerAc = [1 1 1 1];        % spatial streams per AC
guardIntervalUs = [0.8 0.8 0.8 0.8]; % GI in microseconds

baseConfig.totalSlots = 1e5;
baseConfig.arrivalProb = [0.02 0.05 0.01 0.005]; % BK, BE, VI, VO
baseConfig.slotTime = 9e-6; % HE/EHT TXOP limits mapped to 9 us slot model
baseConfig.phyPreambleUs = 100; % HE preamble duration (us)
baseConfig.macHeaderBits = 272; % MAC overhead per PPDU
baseConfig.mcsIndex = mcsPerAc;
baseConfig.nss = nssPerAc;
baseConfig.guardIntervalUs = guardIntervalUs;
baseConfig.acParams = struct( ...
    'name',  {'AC_BK','AC_BE','AC_VI','AC_VO'}, ...
    'aifsn', {7, 3, 2, 2}, ...
    'cwMin', {15, 15, 7, 3}, ...
    'cwMax', {1023, 1023, 15, 7}, ...
    'txopUs', {2528, 2528, 4096, 2080}); % HE TXOP limits in microseconds

for idx = 1:numSweeps
    cfg = baseConfig;
    cfg.numStations = userCounts(idx);
    results = edca_simulation(cfg);

    % Overall (all stations) access success rate per AC.
    acSuccessRate(idx, :) = results.acSuccessRate;

    % Throughput across all ACs for all users.
    acThroughput(idx, :) = results.throughputMbps';
    totalThroughput(idx) = sum(results.throughputMbps);

    % Medium occupancy split between successful PPDUs and collisions.
    busySuccess(idx) = results.busySlotsSuccess / cfg.totalSlots;
    busyCollision(idx) = results.busySlotsCollision / cfg.totalSlots;
end

figure(3); clf;
plot(userCounts, acSuccessRate(:, 1), '-o', ...
     userCounts, acSuccessRate(:, 2), '-s', ...
     userCounts, acSuccessRate(:, 3), '-^', ...
     userCounts, acSuccessRate(:, 4), '-d');
xlabel('Number of stations');
ylabel('Access success rate (successes / attempts)');
title('EDCA Access Success Rate vs. User Count');
legend({'AC_{BK}','AC_{BE}','AC_{VI}','AC_{VO}'}, 'Location', 'southwest', 'Interpreter', 'latex');
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
legend({'Total','AC_{BK}','AC_{BE}','AC_{VI}','AC_{VO}'}, 'Location', 'northwest', 'Interpreter', 'latex');
grid on;

figure(5); clf;
bar(userCounts, [busySuccess busyCollision], 'stacked');
xlabel('Number of stations');
ylabel('Medium time fraction');
title('Medium Occupancy vs. User Count (Success vs. Collision)');
legend({'Success PPDU','Collisions'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
ylim([0 1]);
grid on;
