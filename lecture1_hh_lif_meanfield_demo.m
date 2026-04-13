%% lecture1_hh_lif_meanfield_demo.m
% Classroom MATLAB demo for Lecture 1:
%   1) How Hodgkin-Huxley (HH) reduces to a passive RC membrane near rest,
%      and how adding threshold + reset produces a leaky integrate-and-fire
%      (LIF) caricature.
%   2) How a recurrent noisy LIF network approaches a mean-field closure as
%      N grows.
%   3) How strong synchrony / correlations break that asynchronous
%      mean-field approximation.
%
% The script is self-contained and uses only core MATLAB.
% Tested conceptually for R2016b+ (local functions at end of script).
%
% Suggested classroom flow:
%   - Run section A and discuss HH -> passive membrane -> LIF.
%   - Run section B and discuss law of large numbers / propagation of chaos.
%   - Run section C and discuss why synchrony invalidates the independence
%     assumption behind the simple mean-field closure.

clear; close all; clc;
rng(1);

%% ------------------------------------------------------------------------
% A. HH -> passive membrane -> LIF
% -------------------------------------------------------------------------

hh = defaultHHParams();

% Subthreshold pulse: HH is well approximated by a passive RC membrane if
% we freeze the gates near rest.
I_sub = @(t) 6.0 .* (t >= 10 & t < 30);   % uA/cm^2 (below threshold)
T_sub = 45;                                % ms
DT_HH = 0.01;                              % ms

[t_sub, V_hh_sub]      = simulateHH(I_sub, T_sub, DT_HH, hh);
[~, V_passive_sub]     = simulatePassiveFromHH(I_sub, T_sub, DT_HH, hh);

% Suprathreshold pulse: same HH model, but now compare against an LIF
% abstraction. The spike shape is discarded; we only keep integration,
% threshold crossing, and reset.
I_supra = @(t) 10.0 .* (t >= 10 & t < 13); % short pulse, one spike in HH
T_supra = 35;                               % ms

[t_supra, V_hh_supra]  = simulateHH(I_supra, T_supra, DT_HH, hh);
lif = lifFromHH(hh);
[~, V_lif_supra, lif_spike_times] = simulateLIF(I_supra, T_supra, DT_HH, lif);

hh_spike_times = crossingTimes(t_supra, V_hh_supra, 0.0);   % upward crossing of 0 mV

fprintf('\n==============================================================\n');
fprintf('Lecture 1 demo: HH -> LIF -> mean field -> synchrony failure\n');
fprintf('==============================================================\n');
fprintf('Passive reduction from HH at rest:\n');
fprintf('  V_rest  = %6.2f mV\n', lif.Vrest);
fprintf('  tau_eff = %6.2f ms  (from frozen conductances)\n', lif.tauEff);
fprintf('  R_eff   = %6.2f (mV)/(uA/cm^2)\n', lif.REff);
if ~isempty(hh_spike_times)
    fprintf('First HH spike time under short pulse:   %6.2f ms\n', hh_spike_times(1));
else
    fprintf('HH did not spike under the chosen suprathreshold pulse.\n');
end
if ~isempty(lif_spike_times)
    fprintf('First LIF spike time under short pulse:  %6.2f ms\n', lif_spike_times(1));
else
    fprintf('LIF did not spike under the chosen suprathreshold pulse.\n');
end
fprintf('Interpretation: freeze HH gates near rest -> passive RC membrane;\n');
fprintf('then add threshold + reset -> LIF caricature of spike generation.\n');

figure('Name', 'Lecture 1A: HH to LIF', 'Color', 'w', 'Position', [80 80 1000 720]);

subplot(3,1,1);
plot(t_sub, arrayfun(I_sub, t_sub), 'LineWidth', 1.6);
ylabel('I_{ext}');
title('Subthreshold pulse: full HH vs frozen-gate passive reduction');
box off;

subplot(3,1,2);
plot(t_sub, V_hh_sub, 'LineWidth', 1.8); hold on;
plot(t_sub, V_passive_sub, '--', 'LineWidth', 1.8);
ylabel('V (mV)');
legend({'HH', 'Passive HH linearization'}, 'Location', 'best');
box off;

subplot(3,1,3);
plot(t_supra, V_hh_supra, 'LineWidth', 1.8); hold on;
plot(t_supra, V_lif_supra, '--', 'LineWidth', 1.8);
xlabel('Time (ms)'); ylabel('V (mV)');
legend({'HH', 'LIF'}, 'Location', 'best');
title('Short suprathreshold pulse: spike waveform discarded, threshold/reset retained');
box off;

%% ------------------------------------------------------------------------
% B. Noisy recurrent LIF network and convergence to mean field as N grows
% -------------------------------------------------------------------------

% Dimensionless LIF network for clean population-level numerics.
% V_reset = 0, V_th = 1. Independent noise keeps the activity asynchronous.
netPars = struct();
netPars.tauM     = 20.0;    % ms
netPars.tauRef   = 2.0;     % ms
netPars.tauS     = 5.0;     % ms   (filtered recurrent activity x)
netPars.Vreset   = 0.0;
netPars.Vth      = 1.0;
netPars.muExt    = 1.30;    % mean external drive
netPars.J        = 0.50;    % recurrent gain onto x
netPars.sigmaInd = 0.10;    % independent noise amplitude
netPars.initMode = 'random';
netPars.initSpread = 0.02;

DT_NET = 0.05;              % ms
T_NET  = 500;               % ms
N_list = [100, 500, 2000];

[mf_t, mf_x, mf_r] = simulateMeanField(T_NET, DT_NET, netPars);

asyncOut = cell(size(N_list));
rmse = zeros(size(N_list));
ss_mask = mf_t >= 200;

for k = 1:numel(N_list)
    asyncOut{k} = simulateLIFNetwork(N_list(k), T_NET, DT_NET, netPars, 'async', false, 0);
    rmse(k) = sqrt(mean((asyncOut{k}.x(ss_mask) - mf_x(ss_mask)).^2));
end

figure('Name', 'Lecture 1B: LIF network to mean field', 'Color', 'w', 'Position', [120 120 1100 700]);

subplot(2,1,1);
plot(mf_t, mf_x, 'k', 'LineWidth', 2.0); hold on;
for k = 1:numel(N_list)
    plot(asyncOut{k}.t, asyncOut{k}.x, 'LineWidth', 1.2);
end
xlabel('Time (ms)'); ylabel('x(t)');
legend([{'Mean field'}, arrayfun(@(N) sprintf('Network, N=%d', N), N_list, 'UniformOutput', false)], 'Location', 'best');
title('Asynchronous regime: empirical synaptic activity approaches the mean-field closure');
box off;

subplot(2,1,2);
loglog(N_list, rmse, 'o-', 'LineWidth', 1.8, 'MarkerSize', 8); hold on;
ref_line = rmse(end) * sqrt(N_list(end) ./ N_list);
loglog(N_list, ref_line, '--', 'LineWidth', 1.5);
xlabel('Network size N'); ylabel('RMSE of x(t) vs mean field');
legend({'Measured error', 'Reference slope N^{-1/2}'}, 'Location', 'southwest');
title('Finite-size fluctuations shrink with N (law-of-large-numbers scaling)');
box off;

fprintf('\nAsynchronous LIF network vs mean field:\n');
for k = 1:numel(N_list)
    fprintf('  N = %4d : RMSE(x_N, x_MF) = %.4f\n', N_list(k), rmse(k));
end
fprintf('As N increases, the filtered population activity x_N(t) becomes\n');
fprintf('better approximated by the deterministic mean-field closure.\n');

%% ------------------------------------------------------------------------
% C. Synchrony breaks the simple mean-field closure
% -------------------------------------------------------------------------

syncPars = netPars;
syncPars.muExt    = 1.20;
syncPars.J        = 1.00;
syncPars.sigmaInd = 0.00;   % remove independent noise
syncPars.initMode = 'narrow';
syncPars.initSpread = 0.01; % close-by initial conditions -> phase locking

N_sync = 1000;
nRaster = 100;

syncOut = simulateLIFNetwork(N_sync, T_NET, DT_NET, syncPars, 'sync', true, nRaster);
[mf_sync_t, mf_sync_x, mf_sync_r] = simulateMeanField(T_NET, DT_NET, syncPars);

% Smooth the population rate for visualization.
rate_window_ms = 2.0;
rate_window_pts = max(1, round(rate_window_ms / DT_NET));
rate_kernel = ones(1, rate_window_pts) / rate_window_pts;
sync_rate_smooth = conv(syncOut.popRateHz, rate_kernel, 'same');
mf_rate_smooth   = conv(1000 * mf_sync_r, rate_kernel, 'same');

figure('Name', 'Lecture 1C: Synchrony breaks mean field', 'Color', 'w', 'Position', [150 150 1100 820]);

subplot(3,1,1);
if isempty(syncOut.rasterTimes)
    text(0.5, 0.5, 'No spikes recorded in raster window.', 'HorizontalAlignment', 'center');
    axis off;
else
    plot(syncOut.rasterTimes, syncOut.rasterIds, 'k.', 'MarkerSize', 6);
    xlabel('Time (ms)'); ylabel('Neuron index');
    title(sprintf('Raster (first %d neurons): vertical stripes indicate synchrony', nRaster));
    xlim([0, T_NET]); ylim([0, nRaster + 1]);
    box off;
end

subplot(3,1,2);
plot(syncOut.t, syncOut.x, 'LineWidth', 1.6); hold on;
plot(mf_sync_t, mf_sync_x, '--', 'LineWidth', 1.8);
xlabel('Time (ms)'); ylabel('x(t)');
legend({'Synchronous network', 'Same mean-field closure'}, 'Location', 'best');
title('Correlated bursts create macroscopic oscillations absent from the asynchronous closure');
box off;

subplot(3,1,3);
plot(syncOut.t, sync_rate_smooth, 'LineWidth', 1.5); hold on;
plot(mf_sync_t, mf_rate_smooth, '--', 'LineWidth', 1.8);
xlabel('Time (ms)'); ylabel('Population rate (Hz)');
legend({'Synchronous network', 'Mean-field rate'}, 'Location', 'best');
box off;

fprintf('\nSynchrony counterexample:\n');
fprintf('  Same 1/N coupling normalization does NOT guarantee the simple\n');
fprintf('  asynchronous mean-field closure remains valid.\n');
fprintf('  Once neurons become strongly correlated / synchronized, the\n');
fprintf('  independence assumption is lost and the closure can fail badly.\n');

%% ------------------------------------------------------------------------
% End of main script. Local functions below.
% -------------------------------------------------------------------------

function hh = defaultHHParams()
    hh.C   = 1.0;      % uF/cm^2
    hh.gNa = 120.0;    % mS/cm^2
    hh.gK  = 36.0;     % mS/cm^2
    hh.gL  = 0.3;      % mS/cm^2
    hh.ENa = 50.0;     % mV
    hh.EK  = -77.0;    % mV
    hh.EL  = -54.4;    % mV
    hh.V0  = -65.0;    % mV
end

function [t, V] = simulatePassiveFromHH(Iext, T, dt, hh)
    nSteps = round(T / dt) + 1;
    t = (0:nSteps-1) * dt;
    V = zeros(1, nSteps);

    Vrest = hh.V0;
    m0 = mInf(Vrest);
    h0 = hInf(Vrest);
    n0 = nInf(Vrest);

    gNa0 = hh.gNa * m0^3 * h0;
    gK0  = hh.gK  * n0^4;
    gTot = hh.gL + gNa0 + gK0;
    Eeff = (hh.gL * hh.EL + gNa0 * hh.ENa + gK0 * hh.EK) / gTot;

    V(1) = Vrest;
    for k = 1:nSteps-1
        I = Iext(t(k));
        dV = (I - gTot * (V(k) - Eeff)) / hh.C;
        V(k+1) = V(k) + dt * dV;
    end
end

function lif = lifFromHH(hh)
    Vrest = hh.V0;
    m0 = mInf(Vrest);
    h0 = hInf(Vrest);
    n0 = nInf(Vrest);

    gNa0 = hh.gNa * m0^3 * h0;
    gK0  = hh.gK  * n0^4;
    gTot = hh.gL + gNa0 + gK0;
    Eeff = (hh.gL * hh.EL + gNa0 * hh.ENa + gK0 * hh.EK) / gTot;

    lif = struct();
    lif.tauEff = hh.C / gTot;         % ms
    lif.REff   = 1.0 / gTot;          % (mV)/(uA/cm^2)
    lif.Vrest  = Eeff;                % mV
    lif.Vreset = Eeff;                % mV
    lif.Vth    = -55.0;               % mV, chosen phenomenologically
    lif.tauRef = 4.0;                 % ms, coarse refractory abstraction
end

function [t, V, spikeTimes] = simulateLIF(Iext, T, dt, lif)
    nSteps = round(T / dt) + 1;
    t = (0:nSteps-1) * dt;
    V = zeros(1, nSteps);
    V(1) = lif.Vrest;
    refCount = 0.0;
    spikeTimes = [];

    for k = 1:nSteps-1
        if refCount > 0
            V(k+1) = lif.Vreset;
            refCount = max(0.0, refCount - dt);
            continue;
        end

        I = Iext(t(k));
        dV = (-(V(k) - lif.Vrest) + lif.REff * I) / lif.tauEff;
        V(k+1) = V(k) + dt * dV;

        if V(k+1) >= lif.Vth
            V(k+1) = 20.0;            % draw a spike marker in the trace
            spikeTimes(end+1) = t(k+1); %#ok<AGROW>
            refCount = lif.tauRef;
        end
    end

    % After a spike marker, bring the next sample to reset for visibility.
    spikeIdx = find(V >= 19.0);
    for idx = spikeIdx
        if idx < nSteps
            V(idx+1) = lif.Vreset;
        end
    end
end

function [t, V] = simulateHH(Iext, T, dt, hh)
    nSteps = round(T / dt) + 1;
    t = (0:nSteps-1) * dt;
    V = zeros(1, nSteps);
    m = zeros(1, nSteps);
    h = zeros(1, nSteps);
    n = zeros(1, nSteps);

    V(1) = hh.V0;
    m(1) = mInf(V(1));
    h(1) = hInf(V(1));
    n(1) = nInf(V(1));

    for k = 1:nSteps-1
        I = Iext(t(k));

        INa = hh.gNa * m(k)^3 * h(k) * (V(k) - hh.ENa);
        IK  = hh.gK  * n(k)^4       * (V(k) - hh.EK);
        IL  = hh.gL                 * (V(k) - hh.EL);

        dV = (I - INa - IK - IL) / hh.C;
        dm = alphaM(V(k)) * (1 - m(k)) - betaM(V(k)) * m(k);
        dh = alphaH(V(k)) * (1 - h(k)) - betaH(V(k)) * h(k);
        dn = alphaN(V(k)) * (1 - n(k)) - betaN(V(k)) * n(k);

        V(k+1) = V(k) + dt * dV;
        m(k+1) = m(k) + dt * dm;
        h(k+1) = h(k) + dt * dh;
        n(k+1) = n(k) + dt * dn;
    end
end

function out = simulateLIFNetwork(N, T, dt, pars, mode, storeRaster, nRaster)
    if nargin < 5 || isempty(mode)
        mode = 'async';
    end
    if nargin < 6
        storeRaster = false;
    end
    if nargin < 7
        nRaster = 0;
    end

    nSteps = round(T / dt) + 1;
    t = (0:nSteps-1) * dt;

    switch lower(pars.initMode)
        case 'random'
            V = pars.Vreset + (pars.Vth - pars.Vreset) * rand(N, 1);
        case 'narrow'
            V = pars.Vreset + 0.25 * (pars.Vth - pars.Vreset) + pars.initSpread * randn(N, 1);
        otherwise
            error('Unknown pars.initMode = %s', pars.initMode);
    end

    refractory = zeros(N, 1);
    x = zeros(1, nSteps);
    popRate = zeros(1, nSteps);   % spikes / (neuron * ms)

    if storeRaster
        rasterTimes = [];
        rasterIds = [];
    else
        rasterTimes = [];
        rasterIds = [];
    end

    for k = 1:nSteps-1
        mu = pars.muExt + pars.J * x(k);

        active = refractory <= 0;
        nActive = sum(active);
        if nActive > 0
            drift = dt / pars.tauM * (-V(active) + mu);

            switch lower(mode)
                case 'async'
                    noise = pars.sigmaInd * sqrt(dt / pars.tauM) * randn(nActive, 1);
                case 'sync'
                    noise = zeros(nActive, 1);
                otherwise
                    error('Unknown simulation mode = %s', mode);
            end

            V(active) = V(active) + drift + noise;
        end

        refractory(~active) = max(0.0, refractory(~active) - dt);

        spiking = (V >= pars.Vth) & (refractory <= 0);
        nSpikes = sum(spiking);

        if nSpikes > 0
            if storeRaster && nRaster > 0
                ids = find(spiking(1:min(nRaster, N)));
                if ~isempty(ids)
                    rasterTimes = [rasterTimes; t(k) * ones(numel(ids), 1)]; %#ok<AGROW>
                    rasterIds = [rasterIds; ids(:)]; %#ok<AGROW>
                end
            end

            V(spiking) = pars.Vreset;
            refractory(spiking) = pars.tauRef;
        end

        x(k+1) = x(k) + dt * (-x(k) / pars.tauS) + nSpikes / N;
        popRate(k) = nSpikes / (N * dt);
    end
    popRate(end) = popRate(end-1);

    out = struct();
    out.t = t;
    out.x = x;
    out.popRate = popRate;
    out.popRateHz = 1000 * popRate;
    out.rasterTimes = rasterTimes;
    out.rasterIds = rasterIds;
end

function [t, x, r] = simulateMeanField(T, dt, pars)
    nSteps = round(T / dt) + 1;
    t = (0:nSteps-1) * dt;
    x = zeros(1, nSteps);
    r = zeros(1, nSteps);  % spikes / ms

    for k = 1:nSteps-1
        mu = pars.muExt + pars.J * x(k);
        r(k) = lifTransferDet(mu, pars);
        x(k+1) = x(k) + dt / pars.tauS * (-x(k) + pars.tauS * r(k));
    end
    r(end) = lifTransferDet(pars.muExt + pars.J * x(end), pars);
end

function r = lifTransferDet(mu, pars)
    % Deterministic LIF transfer function (used as a simple classroom-level
    % mean-field closure in the asynchronous regime).
    if mu <= pars.Vth + 1e-12
        r = 0.0;
        return;
    end

    numer = mu - pars.Vreset;
    denom = mu - pars.Vth;

    if numer <= 0 || denom <= 0
        r = 0.0;
        return;
    end

    r = 1.0 / (pars.tauRef + pars.tauM * log(numer / denom));
end

function times = crossingTimes(t, x, threshold)
    idx = find(x(1:end-1) < threshold & x(2:end) >= threshold);
    times = t(idx);
end

function y = alphaN(V)
    y = 0.01 * vtrap(-(V + 55.0), 10.0);
end

function y = betaN(V)
    y = 0.125 * exp(-(V + 65.0) / 80.0);
end

function y = alphaM(V)
    y = 0.1 * vtrap(-(V + 40.0), 10.0);
end

function y = betaM(V)
    y = 4.0 * exp(-(V + 65.0) / 18.0);
end

function y = alphaH(V)
    y = 0.07 * exp(-(V + 65.0) / 20.0);
end

function y = betaH(V)
    y = 1.0 ./ (1.0 + exp(-(V + 35.0) / 10.0));
end

function y = mInf(V)
    y = alphaM(V) ./ (alphaM(V) + betaM(V));
end

function y = hInf(V)
    y = alphaH(V) ./ (alphaH(V) + betaH(V));
end

function y = nInf(V)
    y = alphaN(V) ./ (alphaN(V) + betaN(V));
end

function y = vtrap(x, yscale)
    % Numerically stable version of x / (exp(x / yscale) - 1)
    if abs(x / yscale) < 1e-6
        y = yscale * (1.0 - x / (2.0 * yscale));
    else
        y = x / (exp(x / yscale) - 1.0);
    end
end
