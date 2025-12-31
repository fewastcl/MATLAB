function results = edca_simulation(simConfig)
%EDCA_SIMULATION Discrete-time EDCA simulation for uplink stations to one AP.
%   RESULTS = EDCA_SIMULATION(SIMCONFIG) simulates a single AP serving
%   multiple stations that contend for the medium using the EDCA rules
%   defined by SIMCONFIG. The simulator is intentionally lightweight so it
%   can run inside limited environments but still exposes the main EDCA
%   dynamics: AIFS-based deferrals, random backoff with exponential
%   increase on collision, TXOP durations, and per-AC statistics.
%
%   Required fields of SIMCONFIG:
%     numStations   - number of contending stations.
%     totalSlots    - number of time slots to simulate.
%     arrivalProb   - 1x4 vector with Bernoulli arrival probabilities for
%                     AC_BE, AC_BK, AC_VI, AC_VO in each slot.
%     acParams      - struct array with fields name, aifsn, cwMin, cwMax,
%                     txopSlots that define each access category.
%   Optional fields:
%     slotTime      - slot duration in seconds (default 9e-6).
%     payloadBytes  - payload size for each packet (default 1500).
%
%   The function returns a RESULTS struct containing throughput in Mbps per
%   AC and per station, average delay per AC, total successes, collisions,
%   and a time series of medium usage.

arguments
    simConfig struct
end

% Default parameters
slotTime = getfield_with_default(simConfig, 'slotTime', 9e-6);
payloadBytes = getfield_with_default(simConfig, 'payloadBytes', 1500);

% Validate inputs
numAC = numel(simConfig.acParams);
if numAC ~= 4
    error('Expected four access categories (BE, BK, VI, VO).');
end
if numel(simConfig.arrivalProb) ~= numAC
    error('arrivalProb must be a vector matching acParams.');
end

% Derived constants
payloadBits = payloadBytes * 8;

% Per-station, per-AC state
packetQueues = cell(simConfig.numStations, numAC);
backoffCounters = -1 * ones(simConfig.numStations, numAC);
cwCurrent = repmat([simConfig.acParams.cwMin], simConfig.numStations, 1);
aifsCountdown = zeros(simConfig.numStations, numAC);

% Statistics holders
successCount = zeros(numAC, 1);
collisions = 0;
packetDelays = cell(numAC, 1);
stationSuccess = zeros(simConfig.numStations, numAC);
stationAttempts = zeros(simConfig.numStations, numAC);
acAttempts = zeros(numAC, 1);
mediumBusySlots = 0; % remaining busy slots for the current transmission
previousSlotBusy = false;
mediumUsage = zeros(simConfig.totalSlots, 1); % 0 idle, 1 success, 2 collision

for slotIdx = 1:simConfig.totalSlots
    % Generate arrivals
    for acIdx = 1:numAC
        for sta = 1:simConfig.numStations
            if rand < simConfig.arrivalProb(acIdx)
                packetQueues{sta, acIdx}(end+1) = slotIdx; %#ok<AGROW>
            end
        end
    end

    % If medium is busy, count down and continue
    if mediumBusySlots > 0
        mediumBusySlots = mediumBusySlots - 1;
        previousSlotBusy = true;
        mediumUsage(slotIdx) = 1; % busy due to a previous success
        continue;
    end

    % Medium just became idle after being busy: reset AIFS counters
    if previousSlotBusy
        for acIdx = 1:numAC
            for sta = 1:simConfig.numStations
                if ~isempty(packetQueues{sta, acIdx})
                    aifsCountdown(sta, acIdx) = max(aifsCountdown(sta, acIdx), simConfig.acParams(acIdx).aifsn);
                end
            end
        end
        previousSlotBusy = false;
    end

    % Update AIFS/backoff counters for contenders with queued packets
    for acIdx = 1:numAC
        for sta = 1:simConfig.numStations
            if isempty(packetQueues{sta, acIdx})
                continue;
            end

            if aifsCountdown(sta, acIdx) > 0
                aifsCountdown(sta, acIdx) = aifsCountdown(sta, acIdx) - 1;
                continue;
            end

            if backoffCounters(sta, acIdx) < 0
                backoffCounters(sta, acIdx) = randi(cwCurrent(sta, acIdx) + 1) - 1;
            end

            if backoffCounters(sta, acIdx) > 0
                backoffCounters(sta, acIdx) = backoffCounters(sta, acIdx) - 1;
            end
        end
    end

    % Identify contenders whose backoff hit zero
    contenders = [];
    for acIdx = 1:numAC
        for sta = 1:simConfig.numStations
            if ~isempty(packetQueues{sta, acIdx}) && aifsCountdown(sta, acIdx) == 0 && backoffCounters(sta, acIdx) == 0
                contenders(end+1, :) = [sta, acIdx]; %#ok<AGROW>
            end
        end
    end

    if isempty(contenders)
        mediumUsage(slotIdx) = 0; % idle slot
        continue;
    end

    if size(contenders, 1) > 1
        % Collision
        collisions = collisions + 1;
        mediumUsage(slotIdx) = 2;
        mediumBusySlots = 1; % collision occupies the channel for one slot

        for idx = 1:size(contenders, 1)
            sta = contenders(idx, 1);
            acIdx = contenders(idx, 2);
            stationAttempts(sta, acIdx) = stationAttempts(sta, acIdx) + 1;
            acAttempts(acIdx) = acAttempts(acIdx) + 1;
            cwCurrent(sta, acIdx) = min((cwCurrent(sta, acIdx) + 1) * 2 - 1, simConfig.acParams(acIdx).cwMax);
            backoffCounters(sta, acIdx) = randi(cwCurrent(sta, acIdx) + 1) - 1;
            aifsCountdown(sta, acIdx) = simConfig.acParams(acIdx).aifsn;
        end
        continue;
    end

    % Successful transmission
    sta = contenders(1, 1);
    acIdx = contenders(1, 2);
    queue = packetQueues{sta, acIdx};
    packetArrival = queue(1);
    packetQueues{sta, acIdx}(1) = [];

    stationAttempts(sta, acIdx) = stationAttempts(sta, acIdx) + 1;
    acAttempts(acIdx) = acAttempts(acIdx) + 1;
    successCount(acIdx) = successCount(acIdx) + 1;
    stationSuccess(sta, acIdx) = stationSuccess(sta, acIdx) + 1;
    % Record per-AC packet delay. Appending explicitly avoids any parsing
    % quirks around ``end+1`` on some MATLAB versions.
    delayIdx = numel(packetDelays{acIdx}) + 1;
    packetDelays{acIdx}(delayIdx) = slotIdx - packetArrival;

    cwCurrent(sta, acIdx) = simConfig.acParams(acIdx).cwMin;
    backoffCounters(sta, acIdx) = -1;
    aifsCountdown(sta, acIdx) = simConfig.acParams(acIdx).aifsn;

    txSlots = max(1, simConfig.acParams(acIdx).txopSlots);
    mediumBusySlots = txSlots;
    mediumUsage(slotIdx) = 1;
    previousSlotBusy = true;
end

% Aggregate statistics
throughputMbps = (successCount * payloadBits) / (simConfig.totalSlots * slotTime) / 1e6;
stationThroughput = (stationSuccess * payloadBits) / (simConfig.totalSlots * slotTime) / 1e6;
avgDelay = cellfun(@(d) mean_or_nan(d), packetDelays);
stationSuccessRate = success_divide(stationSuccess, stationAttempts);
acSuccessRate = success_divide(successCount, acAttempts);

results = struct();
results.successCount = successCount;
results.stationSuccess = stationSuccess;
results.stationAttempts = stationAttempts;
results.stationSuccessRate = stationSuccessRate;
results.acAttempts = acAttempts;
results.acSuccessRate = acSuccessRate;
results.throughputMbps = throughputMbps;
results.stationThroughput = stationThroughput;
results.avgDelaySlots = avgDelay;
results.collisions = collisions;
results.mediumUsage = mediumUsage;
results.config = simConfig;
end

function val = getfield_with_default(s, fieldName, defaultValue)
if isfield(s, fieldName)
    val = s.(fieldName);
else
    val = defaultValue;
end
end

function out = mean_or_nan(vec)
if isempty(vec)
    out = NaN;
else
    out = mean(vec);
end
end

function rate = success_divide(successes, attempts)
rate = successes ./ attempts;
rate(attempts == 0) = NaN;
end
