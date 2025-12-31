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
%     mcsIndex      - scalar or 1x4 vector of UHR-MCS indices (0-based) per AC.
%     nss           - scalar or 1x4 vector of spatial streams per AC (default 1).
%     guardIntervalUs - scalar or 1x4 vector GI in microseconds (0.8/1.6/3.2).
%     phyRateMbps   - 1x4 vector override of PHY data rates per AC.
%     phyPreambleUs - PHY preamble duration in microseconds (default 100 us).
%     macHeaderBits - MAC overhead bits per PPDU (default 272 bits ~34 bytes).
%     collisionDurationUs - minimum medium busy time for a collision (default preamble).
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
phyPreambleUs = getfield_with_default(simConfig, 'phyPreambleUs', 100);
macHeaderBits = getfield_with_default(simConfig, 'macHeaderBits', 272);
collisionDurationUs = getfield_with_default(simConfig, 'collisionDurationUs', phyPreambleUs);
mcsIndex = expand_to_ac(getfield_with_default(simConfig, 'mcsIndex', [3 3 7 9]), numAC);
nss = expand_to_ac(getfield_with_default(simConfig, 'nss', 1), numAC);
guardIntervalUs = expand_to_ac(getfield_with_default(simConfig, 'guardIntervalUs', 0.8), numAC);

if isfield(simConfig, 'phyRateMbps')
    phyRateMbps = simConfig.phyRateMbps;
else
    phyRateMbps = zeros(1, numAC);
    for acIdx = 1:numAC
        phyRateMbps(acIdx) = uhr_phy_rate(mcsIndex(acIdx), nss(acIdx), guardIntervalUs(acIdx));
    end
end

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
collisionBusySlots = max(1, ceil((collisionDurationUs * 1e-6) / slotTime));
txopLimitSlots = ones(1, numAC);
for acIdx = 1:numAC
    if isfield(simConfig.acParams(acIdx), 'txopUs')
        txopLimitSlots(acIdx) = max(1, ceil(simConfig.acParams(acIdx).txopUs * 1e-6 / slotTime));
    elseif isfield(simConfig.acParams(acIdx), 'txopSlots')
        txopLimitSlots(acIdx) = max(1, simConfig.acParams(acIdx).txopSlots);
    else
        txopLimitSlots(acIdx) = 1;
    end
end

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
busyCause = 0; % 1 = success, 2 = collision
busySlotsSuccess = 0;
busySlotsCollision = 0;
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
        if busyCause == 1
            busySlotsSuccess = busySlotsSuccess + 1;
        elseif busyCause == 2
            busySlotsCollision = busySlotsCollision + 1;
        end
        mediumBusySlots = mediumBusySlots - 1;
        previousSlotBusy = true;
        mediumUsage(slotIdx) = busyCause; % track current busy reason
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
        % Collision duration: attempted PPDU (capped by TXOP), at least preamble.
        collidedTx = zeros(size(contenders, 1), 1);
        for cIdx = 1:size(contenders, 1)
            acIdx = contenders(cIdx, 2);
            ppduSlots = compute_ppdu_slots(payloadBits, macHeaderBits, phyPreambleUs, phyRateMbps(acIdx), slotTime);
            collidedTx(cIdx) = min(ppduSlots, txopLimitSlots(acIdx));
        end
        mediumBusySlots = max(collisionBusySlots, max(collidedTx));
        busyCause = 2;

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

    % Compute PPDU duration based on PHY rate and overhead, capped by TXOP.
    ppduSlots = compute_ppdu_slots(payloadBits, macHeaderBits, phyPreambleUs, phyRateMbps(acIdx), slotTime);
    txSlots = min(ppduSlots, txopLimitSlots(acIdx));
    mediumBusySlots = txSlots;
    mediumUsage(slotIdx) = 1;
    busyCause = 1;
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
results.busySlotsSuccess = busySlotsSuccess;
results.busySlotsCollision = busySlotsCollision;
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

function rateMbps = uhr_phy_rate(mcsIndex, nss, guardIntervalUs)
%UHR_PHY_RATE Return data rate in Mbps for given MCS, NSS, GI (us) using 242-tone RU table.
ndbpsTable = [117 351 468 702 936 1053 1170 1404 1755 1950 2340 59 234 468 702 936 1404];
if mcsIndex < 0 || mcsIndex >= numel(ndbpsTable)
    error('mcsIndex must be between 0 and %d.', numel(ndbpsTable)-1);
end
ndbps = ndbpsTable(mcsIndex + 1) * max(1, nss);
tsym = 12.8e-6 + guardIntervalUs * 1e-6;
rateMbps = (ndbps / tsym) / 1e6;
end

function vec = expand_to_ac(val, numAC)
if isscalar(val)
    vec = repmat(val, 1, numAC);
else
    vec = val;
end
end

function slots = compute_ppdu_slots(payloadBits, macHeaderBits, preambleUs, phyRateMbps, slotTime)
dataBits = payloadBits + macHeaderBits;
txSeconds = dataBits / (phyRateMbps * 1e6);
totalSeconds = txSeconds + preambleUs * 1e-6;
slots = max(1, ceil(totalSeconds / slotTime));
end
