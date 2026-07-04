%% PARAMETERS

mod_depth = 0.1;        % Modulation Depth
fc = 950e6;             % Carrier Frequency
fs = 1e6;               % Sampling Frequency
sps = 12;               % Samples per Symbol ( for upsampling )
M = 4;                  % Modulation Index

Ndata = 50000;          % Length of data symbols

rolloff = 0.35;         % Roll-off factor of RRC Filter
span = 8;               % Span of RRC Filter

% Local ZC Sequence

Nsync = 1021;           % Length of Sync Data
u = 41;                 % Factor in Zadoff Chu Sequence

% Nsync & u must be co-prime

n = (0:Nsync-1).';

zcSync = exp(-1j*pi*u*n.*(n+1)/Nsync);

% Pilot Data

Npilot = 7500;          % Length of Pilot Symbols
rng(123);               % Using rng(seed) to use the exact same pilot data on both TX & RX ends.

pilotBits = randi([0 1], Npilot, 1);
pilotSym  = pskmod(pilotBits,2,0,'gray');

% REGENERATE DATA

rng(12345);            % Using rng(seed) to use the exact same data on both TX & RX ends.

dataBits = randi([0 M-1], Ndata, 1);
dataSym  = mod_depth * qammod(dataBits, M, 'UnitAveragePower', true);

    
frameLengthSym = (Nsync + Npilot + Ndata);

frameSamples = ...
    (frameLengthSym-1)*sps + span*sps + 1;          % Addition of RRC Filter length

% %% Frame

txSym = [zcSync; pilotSym; dataSym];

rrcFilter = rcosdesign(rolloff,span,sps,'sqrt');

txSignal = upfirdn(txSym,rrcFilter,sps,1);


%% RX

rx = sdrrx('Pluto');

rx.CenterFrequency = fc;

rx.BasebandSampleRate = fs;

rx.GainSource = 'Manual';

rx.Gain = 0;            % This is fixed after lots of trials, and gain saturation observations

rx.OutputDataType = 'double';

rx.SamplesPerFrame = frameSamples*2;

%% Transmit and Release

rxRaw = rx();
pause(2);

release(rx);


%% RRC & Coarse CFO

rrc = rcosdesign(rolloff,span,sps,'sqrt');

y = rxRaw.^M;

phaseDiff = angle( ...
    y(2:end).*conj(y(1:end-1)));

cfoEst = ...
angle(mean( ...
y(2:end).*conj(y(1:end-1)))) ...
* fs/(2*M*pi);

fprintf('Coarse CFO = %.2f Hz\n',cfoEst);

fprintf('Norm RX = %.3f\n', ...
    mean(abs(rxRaw).^2));

n = (0:length(rxRaw)-1).';

rxRaw = rxRaw .* ...
    exp(-1j*2*pi*cfoEst*n/fs);

%% MATCHED FILTER & FRAME SYNC.

rxMF = conv(rxRaw,rrc,'same');

syncWave = upfirdn(zcSync,rrc,sps,1);

syncMF = conv(syncWave,rrc,'same');

[corrOut,lags] = xcorr(rxMF,syncMF);

[pks,locs] = findpeaks(abs(corrOut));

[pks,I] = sort(pks,'descend');

peak1 = pks(1);
peak2 = pks(2);

lag1 = lags(locs(I(1)));
lag2 = lags(locs(I(2)));

if peak1/peak2 < 1.1
    lag = min(lag1,lag2);
else
    lag = lag1;
end


%%
figure;
plot(abs(corrOut))
title('Correlation Peaks')          
% If there are too many peaks, 
% either data is not received properly or, 
% sync data must be increased in length

%% PACKET EXTRACTION

frameLengthSym = Nsync + Npilot + Ndata;

waveLength = ...
    (frameLengthSym-1)*sps + length(rrc);

startIndex = lag + 1;

endIndex = startIndex + waveLength - 1;

if endIndex <= length(rxMF)

    rxFrame = rxMF(startIndex:endIndex);

else

    error('Frame exceeds capture');

end



%% DOWNSAMPLE
half = sps/2;
delay = sps*span/2;         % Delay due to RRC Filtering


%% Obtaining the best offset using the Pilot Metric using Correlation
% ( Better the correlation, better offset )
metric = zeros(sps,1);
for offset = 1:sps

    rxSymTest = rxFrame(delay+offset:sps:end);

    if length(rxSymTest) < Nsync+Npilot
        continue;
    end

    rxPilot = rxSymTest(Nsync+1:Nsync+Npilot);

    metric(offset) = ...
        abs(sum(rxPilot .* conj(pilotSym)));

end


[~,bestOffset] = max(metric);

fprintf('Best offset = %d\n',bestOffset);

rxSym = rxFrame(delay+bestOffset:sps:end);

L = min(length(rxSym),frameLengthSym);

rxSym = rxSym(1:L);


%% CFO CORRECTION

rxSymCFO = rxSym;
for iter = 1:100

    if length(rxSymCFO) < Nsync + Npilot
        error('Pilot not fully received');
    end

    rxPilot = rxSymCFO(Nsync+1:Npilot + Nsync);

    err = rxPilot .* conj(pilotSym);

    phaseErr = unwrap(angle(err));

    nPilot = (0:length(phaseErr)-1).';

    p = polyfit(nPilot,phaseErr,1);

    slope = p(1);

    fprintf('Iter %d : slope = %.4e\n', ...
        iter,slope);

    if abs(slope) < 1e-7        % Can be adjusted with respect to data length
        break;
    end

    nAll = (0:length(rxSymCFO)-1).';

    rxSymCFO = ...
        rxSymCFO .* exp(-1j*slope*nAll);

end
pause(1);


%% RESIDUAL PHASE
if length(rxSymCFO) < Nsync + Npilot
    error('Pilot not fully received');
end
rxPilot = rxSymCFO(Nsync+1:Npilot + Nsync);

err = rxPilot .* conj(pilotSym);

phaseErr = unwrap(angle(err));

p = polyfit(nPilot,phaseErr,1);

rxPilot = rxSymCFO(Nsync+1:Nsync+Npilot);

phaseOffsetEst = angle( ...
    sum(rxPilot .* conj(pilotSym)));

rxSymCFO = rxSymCFO .* ...
    exp(-1j*phaseOffsetEst);

rxPilot = rxSymCFO(Nsync+1:Nsync+Npilot);

fprintf('Residual Slope = %.4e\n',p(1));
fprintf('Phase error = %.4e\n',mean(phaseErr));
fprintf('Phase offset est = %.4e\n',phaseOffsetEst);
pause(1);



%% CHANNEL ESTIMATION

h_est = sum(rxPilot .* conj(pilotSym )) ...
    / sum(abs(pilotSym).^2);

fprintf('|h_est| = %.3f\n',abs(h_est));

fprintf('Phase(h_est) = %.2f deg\n', ...
    angle(h_est)*180/pi);

pause(1);

figure;
plot(real(rxSymCFO),imag(rxSymCFO),'.');
axis equal;
grid on;

%% EQUALIZATION

eqSym = rxSymCFO ./ h_est;
% mean(eqSym)
rxData = eqSym(Nsync+Npilot+1:end);

% Removing DC offset from only data symbols on 2 iterations
rxData = rxData - mean(rxData);

rxData = rxData - mean(rxData);

% Trim to same length as bit vectors
Nsym = min(length(rxData), Ndata);

dataBits   = dataBits(1:Nsym); 

rxDataBits = qamdemod(rxData, M, 'UnitAveragePower', true);
rxDataBits = rxDataBits(1:Nsym);

%%

figure;
plot(real(rxData),imag(rxData),'.');
axis equal;
grid on;

%%

% Build reference constellation
refConst = mod_depth * qammod((0:M-1).', M, 'UnitAveragePower', true);

% Euclidean distances
distMat = abs(rxData - refConst.');
[minDist, nearestIdx] = min(distMat, [], 2);

rxDecoded = nearestIdx - 1;

% Decision radius with respect to mod_depth
a = mod_depth;
d_min = sqrt(2) * a;
r_tangent = 0.6*(d_min/2);

% Decision radius with respect to sigma
noise = rxData - refConst(nearestIdx);
sigma = sqrt(mean(abs(noise).^2));
rho = 0.35;
r_rho = rho*sigma;

% Choose confidence radius
r_decision = min(r_tangent,r_rho);

% Apply same radius to every constellation point
validSym = minDist < r_decision;

% SER
ser = mean(rxDecoded(validSym) ~= dataBits(validSym));

fprintf('Accepted(Radial Gating) = %d / %d\n', sum(validSym), Nsym);
fprintf('SER(Radial Gating) = %.6f\n', ser);
%%
figure; hold on;

% Masks based on transmitted symbols
sym0 = (dataBits == 0);
sym1 = (dataBits == 1);
sym2 = (dataBits == 2);
sym3 = (dataBits == 3);

% Plot received symbols
plot(real(rxData(sym0)), imag(rxData(sym0)), '.', ...
    'Color', [1 0 0], 'MarkerSize', 5);      % Red

plot(real(rxData(sym1)), imag(rxData(sym1)), '.', ...
    'Color', [0 0.6 1], 'MarkerSize', 5);    % Blue

plot(real(rxData(sym2)), imag(rxData(sym2)), '.', ...
    'Color', [0 0.7 0], 'MarkerSize', 5);    % Green

plot(real(rxData(sym3)), imag(rxData(sym3)), '.', ...
    'Color', [1 0.5 0], 'MarkerSize', 5);    % Orange

% Overlay ideal constellation points
plot(real(refConst(1)), imag(refConst(1)), 'r+', ...
    'MarkerSize', 12, 'LineWidth', 2);

plot(real(refConst(2)), imag(refConst(2)), 'b+', ...
    'MarkerSize', 12, 'LineWidth', 2);

plot(real(refConst(3)), imag(refConst(3)), 'g+', ...
    'MarkerSize', 12, 'LineWidth', 2);

plot(real(refConst(4)), imag(refConst(4)), ...
    '+', 'Color', [1 0.5 0], ...
    'MarkerSize', 12, 'LineWidth', 2);

% Draw decision circles (optional)
theta = linspace(0,2*pi,100);

for k = 1:4
    cx = real(refConst(k));
    cy = imag(refConst(k));

    plot(cx + r_decision*cos(theta), ...
        cy + r_decision*sin(theta), ...
        'w', 'LineWidth', 0.8);
end

axis equal;
grid on;
xlabel('I');
ylabel('Q');
title('4-QAM Constellation');

legend('Symbol 0','Symbol 1','Symbol 2','Symbol 3', ...
    'Ideal 0','Ideal 1','Ideal 2','Ideal 3', ...
    'Location','best');


%% Rectangular decision regions

% Build reference constellation
refConst = mod_depth * qammod((0:M-1).', M, 'UnitAveragePower', true);

% Find nearest constellation point
distMat = abs(rxData - refConst.');
[~, nearestIdx] = min(distMat, [], 2);

rxDecoded = nearestIdx - 1;

% Rectangle half-width

a = mod_depth;

% (or adaptive)
noise = rxData - refConst(nearestIdx);
sigma = sqrt(mean(abs(noise).^2));
rectHalf = min(r_tangent, rho*sigma);

% Compute I/Q errors

Ierr = abs(real(rxData) - real(refConst(nearestIdx)));
Qerr = abs(imag(rxData) - imag(refConst(nearestIdx)));

% Accept only points inside rectangle

validSym = (Ierr < rectHalf) & (Qerr < rectHalf);

% SER

ser = mean(rxDecoded(validSym) ~= dataBits(validSym));

fprintf('Accepted(Rectangular Gating) = %d / %d\n', sum(validSym), Nsym);
fprintf('SER(Rectangular Gating) = %.6f\n', ser);
%%
figure; hold on;

% Masks based on transmitted symbol
sym0 = (dataBits == 0);
sym1 = (dataBits == 1);
sym2 = (dataBits == 2);
sym3 = (dataBits == 3);

% Plot received symbols

plot(real(rxData(sym0)), imag(rxData(sym0)), '.', ...
    'Color',[0 0.4470 0.7410], 'MarkerSize',4);     % Blue

plot(real(rxData(sym1)), imag(rxData(sym1)), '.', ...
    'Color',[0.8500 0.3250 0.0980], 'MarkerSize',4); % Orange

plot(real(rxData(sym2)), imag(rxData(sym2)), '.', ...
    'Color',[0.4660 0.6740 0.1880], 'MarkerSize',4); % Green

plot(real(rxData(sym3)), imag(rxData(sym3)), '.', ...
    'Color',[0.4940 0.1840 0.5560], 'MarkerSize',4); % Purple

% Ideal constellation points

plot(real(refConst(1)), imag(refConst(1)), ...
    'b+', 'MarkerSize',12,'LineWidth',2);

plot(real(refConst(2)), imag(refConst(2)), ...
    'r+', 'MarkerSize',12,'LineWidth',2);

plot(real(refConst(3)), imag(refConst(3)), ...
    'g+', 'MarkerSize',12,'LineWidth',2);

plot(real(refConst(4)), imag(refConst(4)), ...
    'm+', 'MarkerSize',12,'LineWidth',2);

% Draw acceptance rectangles

for k = 1:4

    cx = real(refConst(k));
    cy = imag(refConst(k));

    rectangle( ...
        'Position', [ ...
            cx-rectHalf,...
            cy-rectHalf,...
            2*rectHalf,...
            2*rectHalf], ...
        'EdgeColor','w', ...
        'LineWidth',1);

end

axis equal;
grid on;
xlabel('I');
ylabel('Q');

title('Received 4-QAM Constellation with Rectangular Decision Regions');

legend( ...
    'Symbol 0 RX', ...
    'Symbol 1 RX', ...
    'Symbol 2 RX', ...
    'Symbol 3 RX', ...
    'Ideal 0', ...
    'Ideal 1', ...
    'Ideal 2', ...
    'Ideal 3', ...
    'Location','best');

%% Reference constellation
refConst = mod_depth * qammod((0:M-1).',M,'UnitAveragePower',true);

% Nearest symbol (for decoding only)

distMat = abs(rxData-refConst.');
[~,nearestIdx] = min(distMat,[],2);

rxDecoded = nearestIdx-1;

% Angular threshold

theta_th = deg2rad(12);      % ±12 degrees

% Radius threshold

r_min = min(r_tangent,r_rho)*7.3;
r_max = 5.55 * mod_depth;

% Phase of received symbols

rxPhase = angle(rxData);

% Ideal phase of nearest constellation point

refPhase = angle(refConst(nearestIdx));

% Wrapped phase error

phaseErr = angle(exp(1j*(rxPhase-refPhase)));

% Magnitude

mag = abs(rxData);

% Accept symbols

validSym = ...
    (abs(phaseErr) < theta_th) & ...
    (mag > r_min);

% SER

ser = mean(rxDecoded(validSym) ~= dataBits(validSym));

fprintf('Accepted (Angular Gating) = %d / %d\n',sum(validSym),Nsym);
fprintf('SER (Angular Gating) = %.6f\n',ser);

%%
figure; hold on;

% Symbol masks
sym0 = (dataBits == 0);
sym1 = (dataBits == 1);
sym2 = (dataBits == 2);
sym3 = (dataBits == 3);

% Received symbols

plot(real(rxData(sym0)), imag(rxData(sym0)), '.', ...
    'Color',[0 0.4470 0.7410], 'MarkerSize',4);

plot(real(rxData(sym1)), imag(rxData(sym1)), '.', ...
    'Color',[0.8500 0.3250 0.0980], 'MarkerSize',4);

plot(real(rxData(sym2)), imag(rxData(sym2)), '.', ...
    'Color',[0.4660 0.6740 0.1880], 'MarkerSize',4);

plot(real(rxData(sym3)), imag(rxData(sym3)), '.', ...
    'Color',[0.4940 0.1840 0.5560], 'MarkerSize',4);

% Ideal constellation

plot(real(refConst(1)), imag(refConst(1)), ...
    'b+','LineWidth',2,'MarkerSize',12);

plot(real(refConst(2)), imag(refConst(2)), ...
    'r+','LineWidth',2,'MarkerSize',12);

plot(real(refConst(3)), imag(refConst(3)), ...
    'g+','LineWidth',2,'MarkerSize',12);

plot(real(refConst(4)), imag(refConst(4)), ...
    'm+','LineWidth',2,'MarkerSize',12);

% Inner and outer circles

theta = linspace(0,2*pi,500);

plot(r_min*cos(theta), ...
     r_min*sin(theta), ...
     'w-','LineWidth',1.5);

plot(r_max*cos(theta), ...
     r_max*sin(theta), ...
     'w-','LineWidth',1.5);

% Angular sector boundaries

for k = 1:4

    phi = angle(refConst(k));

    % Lower boundary
    plot( ...
        [r_min*cos(phi-theta_th), ...
         r_max*cos(phi-theta_th)], ...
        [r_min*sin(phi-theta_th), ...
         r_max*sin(phi-theta_th)], ...
        'r','LineWidth',1);

    % Upper boundary
    plot( ...
        [r_min*cos(phi+theta_th), ...
         r_max*cos(phi+theta_th)], ...
        [r_min*sin(phi+theta_th), ...
         r_max*sin(phi+theta_th)], ...
        'r','LineWidth',1);

    % Inner arc
    t = linspace(phi-theta_th,phi+theta_th,100);

    plot( ...
        r_min*cos(t), ...
        r_min*sin(t), ...
        'w','LineWidth',1);

    % Outer arc
    plot( ...
        r_max*cos(t), ...
        r_max*sin(t), ...
        'w','LineWidth',1);

end

% Draw acceptance rectangles

for k = 1:4

    cx = real(refConst(k));
    cy = imag(refConst(k));

    rectangle( ...
        'Position', [ ...
        cx-rectHalf,...
        cy-rectHalf,...
        2*rectHalf,...
        2*rectHalf], ...
        'EdgeColor','y', ...
        'LineWidth',1);

end

% Draw decision circles (optional)
theta = linspace(0,2*pi,100);

for k = 1:4
    cx = real(refConst(k));
    cy = imag(refConst(k));

    plot(cx + r_decision*cos(theta), ...
        cy + r_decision*sin(theta), ...
        'g', 'LineWidth', 0.8);
end

axis equal;
grid on;

xlabel('I');
ylabel('Q');

title('4-QAM Cartesian-Angular-Radial Decision Regions');

legend( ...
    'Symbol 0', ...
    'Symbol 1', ...
    'Symbol 2', ...
    'Symbol 3', ...
    'Ideal 0', ...
    'Ideal 1', ...
    'Ideal 2', ...
    'Ideal 3', ...
    'Location','best');
%%

txBits = de2bi(dataBits(validSym),2,'left-msb');
rxBits = de2bi(rxDecoded(validSym),2,'left-msb');

QBER = sum(txBits(:)~=rxBits(:))/numel(txBits);

fprintf("QBER = %.6f\n",QBER);