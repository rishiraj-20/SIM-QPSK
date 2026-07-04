%% DSP Function
function out = receiverDSP( ...
    rxRaw,...
    M,...
    fs,...
    rolloff,...
    span,...
    sps,...
    zcSync,...
    pilotSym,...
    Nsync,...
    Npilot,...
    Ndata)

%% ===========================
%% RRC & Coarse CFO
%% ===========================

rrc = rcosdesign(rolloff,span,sps,'sqrt');

y = rxRaw.^M;

out.cfoEst = angle(mean( ...
    y(2:end).*conj(y(1:end-1)))) ...
    * fs/(2*M*pi);

fprintf('Coarse CFO = %.2f Hz\n',out.cfoEst);

fprintf('Norm RX = %.3f\n', ...
    mean(abs(rxRaw).^2));

n = (0:length(rxRaw)-1).';

rxRaw = rxRaw .* ...
    exp(-1j*2*pi*out.cfoEst*n/fs);

%% ===========================
%% MATCHED FILTER & FRAME SYNC
%% ===========================

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

if peak1/peak2 < 1.3
    lag = min(lag1,lag2);
else
    lag = lag1;
end

out.peakRatio = peak1/peak2;
out.lag = lag;

%% ===========================
%% PACKET EXTRACTION
%% ===========================

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

%% ===========================
%% DOWNSAMPLE
%% ===========================

delay = sps*span/2;

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

out.bestOffset = bestOffset;
out.metric = metric;

rxSym = rxFrame(delay+bestOffset:sps:end);

L = min(length(rxSym),frameLengthSym);

rxSym = rxSym(1:L);

%% ===========================
%% CFO CORRECTION
%% ===========================

rxSymCFO = rxSym;

for iter = 1:100

    if length(rxSymCFO) < Nsync+Npilot
        error('Pilot not fully received');
    end

    rxPilot = rxSymCFO(Nsync+1:Npilot+Nsync);

    err = rxPilot .* conj(pilotSym);

    phaseErr = unwrap(angle(err));

    nPilot = (0:length(phaseErr)-1).';

    p = polyfit(nPilot,phaseErr,1);

    slope = p(1);

    fprintf('Iter %d : slope = %.4e\n',iter,slope);

    if abs(slope) < 1e-16
        break;
    end

    nAll = (0:length(rxSymCFO)-1).';

    rxSymCFO = ...
        rxSymCFO .* exp(-1j*slope*nAll);

end

%% ===========================
%% RESIDUAL PHASE
%% ===========================

rxPilot = rxSymCFO(Nsync+1:Nsync+Npilot);

err = rxPilot .* conj(pilotSym);

phaseErr = unwrap(angle(err));

p = polyfit(nPilot,phaseErr,1);

phaseOffsetEst = angle( ...
    sum(rxPilot .* conj(pilotSym)));

rxSymCFO = rxSymCFO .* ...
    exp(-1j*phaseOffsetEst);

rxPilot = rxSymCFO(Nsync+1:Nsync+Npilot);

fprintf('Residual Slope = %.4e\n',p(1));
fprintf('Phase error = %.4e\n',mean(phaseErr));
fprintf('Phase offset est = %.4e\n',phaseOffsetEst);

out.phaseSlope = p(1);
out.phaseOffset = phaseOffsetEst;

%% ===========================
%% CHANNEL ESTIMATION
%% ===========================

h_est = ...
    sum(rxPilot .* conj(pilotSym)) ...
    / sum(abs(pilotSym).^2);

fprintf('|h_est| = %.3f\n',abs(h_est));

fprintf('Phase(h_est) = %.2f deg\n', ...
    angle(h_est)*180/pi);

out.h_est = h_est;

%% ===========================
%% EQUALIZATION
%% ===========================

eqSym = rxSymCFO ./ h_est;

rxData = eqSym(Nsync+Npilot+1:end);

out.meanFirstHalf = mean(rxData(1:5000));
out.meanSecondHalf = mean(rxData(5001:end));

rxData = rxData - mean(rxData);
rxData = rxData - mean(rxData);

out.meanFinal = mean(rxData);

%% ===========================
%% OUTPUTS
%% ===========================

out.rrc = rrc;
out.rxMF = rxMF;
out.rxFrame = rxFrame;

out.rxSym = rxSym;
out.rxSymCFO = rxSymCFO;

out.eqSym = eqSym;

out.rxPilot = rxPilot;

out.rxData = rxData;

end

%% Mutual Information Function
function I = MutualInfo(tx,rx,M)

N = length(tx);

joint = zeros(M,M);

for k = 1:N
    joint(tx(k)+1,rx(k)+1) = ...
        joint(tx(k)+1,rx(k)+1) + 1;
end

joint = joint/N;

Px = sum(joint,2);
Py = sum(joint,1);

I = 0;

for i = 1:M
    for j = 1:M

        if joint(i,j)>0

            I = I + ...
                joint(i,j)* ...
                log2(joint(i,j)/(Px(i)*Py(j)));

        end

    end
end

end
%% PARAMETERS

rhoVec = 0.05:0.05:5;       % Vector for rho sweep
Nrho = length(rhoVec);

thetaVec = 1:1:45;          % Vector for theta sweep
Nt = length(thetaVec);

% All required Parameters for Comparative Analysis

% Accepted number of symbols
accepted_rad_B  = zeros(Nrho,1);
accepted_rect_B = zeros(Nrho,1);
accepted_ang_B  = zeros(Nrho,1);

% Psift 
Pacc_rad  = zeros(Nrho,1);
Pacc_rect = zeros(Nrho,1);
Pacc_ang  = zeros(Nrho,1);

% SER of different detectors
ser_rad_B  = nan(Nrho,1);
ser_rect_B = nan(Nrho,1);
ser_ang_B  = nan(Nrho,1);

% Secret Fraction
SecretFraction_rad  = nan(Nrho,1);
SecretFraction_rect = nan(Nrho,1);
SecretFraction_ang  = nan(Nrho,1);

% ESKR 
ESKR_rad_kbps  = nan(Nrho,1);
ESKR_rect_kbps = nan(Nrho,1);
ESKR_ang_kbps  = nan(Nrho,1);

% BER of Eve
BER_eve_rad  = nan(Nrho,1);
BER_eve_rect = nan(Nrho,1);
BER_eve_ang  = nan(Nrho,1);

% SER of Eve
SER_eve_rad  = nan(Nrho,1);
SER_eve_rect = nan(Nrho,1);
SER_eve_ang  = nan(Nrho,1);

% QBER
QBER_rad  = nan(Nrho,1);
QBER_rect = nan(Nrho,1);
QBER_ang  = nan(Nrho,1);

mod_depth = 0.1;        % Modulation Depth

fc = 950e6;             % Carrier Frequency
fs = 1e6;               % Sampling Frequency
sps = 12;               % Samples per Symbol (for upsampling)
M = 4;                  % Modulation Index

Ndata = 5e4;            % Length of Data Symbol

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

%% Frame

txSym = [zcSync; pilotSym; dataSym];

rrcFilter = rcosdesign(rolloff,span,sps,'sqrt');

txSignal = upfirdn(txSym,rrcFilter,sps,1);


%% RX Bob and Eve

rxBob = sdrrx('Pluto');
rxBob.RadioID = "usb:0";        % USB Address for Bob RX

rxBob.CenterFrequency = fc;
rxBob.BasebandSampleRate = fs;
rxBob.GainSource = 'Manual';
rxBob.Gain = 0;         % This is fixed after lots of trials, and gain saturation observations
rxBob.OutputDataType = 'double';
rxBob.SamplesPerFrame = frameSamples*2;


rxEve = sdrrx('Pluto');
rxEve.RadioID = "usb:1";        % USB Address for Eve RX

rxEve.CenterFrequency = fc;
rxEve.BasebandSampleRate = fs;
rxEve.GainSource = 'Manual';
rxEve.Gain = 0;         % This is fixed after lots of trials, and gain saturation observations
rxEve.OutputDataType = 'double';
rxEve.SamplesPerFrame = frameSamples*2;

%%

disp("Receiving...");

rxRaw_B = rxBob();
rxRaw_E = rxEve();

disp("Reception Complete");

release(rxBob);
release(rxEve);

Bob = receiverDSP( ...
    rxRaw_B,...
    M,...
    fs,...
    rolloff,...
    span,...
    sps,...
    zcSync,...
    pilotSym,...
    Nsync,...
    Npilot,...
    Ndata);

Eve = receiverDSP( ...
    rxRaw_E,...
    M,...
    fs,...
    rolloff,...
    span,...
    sps,...
    zcSync,...
    pilotSym,...
    Nsync,...
    Npilot,...
    Ndata);

%%
rxData_B = Bob.rxData;
rxData_B = rxData_B(:);
dataBits_B = dataBits;
h_est = Bob.h_est;
rxData_E = Eve.rxData;


refConst = mod_depth * qammod(0:M-1,M,'UnitAveragePower',true);
refConst = refConst(:);
dist = abs(rxData_B - refConst.');

[minDist, nearestIdx] = min(dist, [], 2);

Ierr = abs(real(rxData_B) - real(refConst(nearestIdx)));
Qerr = abs(imag(rxData_B) - imag(refConst(nearestIdx)));

phaseErr1 = angle(rxData_B .* conj(refConst(nearestIdx)));
mag = abs(rxData_B);

rxDecoded_B = nearestIdx-1;

rxDecoded_E = qamdemod( ...
    rxData_E/mod_depth,...
    M,...
    'UnitAveragePower',true);

noise = rxData_B - refConst(nearestIdx);

sigma_B = sqrt(mean(abs(noise).^2));

fprintf('Noise sigma = %.4f\n', sigma_B);

Nsym = Ndata;

%% 
for ii=1:Nrho

    rho = rhoVec(ii);

    %% =========================
    %% BOB
    %% =========================

    r_decision = rho*sigma_B;

    %% Radial

    valid_rad = minDist < r_decision;

    accepted_rad_B(ii)=sum(valid_rad);

    if any(valid_rad)

        tx_B = dataBits_B(valid_rad);
        rx_B = rxDecoded_B(valid_rad);

        I_AB = MutualInfo(tx_B,rx_B,M);

        tx_E = dataBits_B(valid_rad);
        rx_E = rxDecoded_E(valid_rad);

        I_AE = MutualInfo(tx_E,rx_E,M);

        ser_rad_B(ii) = mean(tx_B~=rx_B);
        SER_eve_rad(ii) = mean(tx_E~=rx_E);

        txBitsB = de2bi(tx_B,2,'left-msb');
        rxBitsB = de2bi(rx_B,2,'left-msb');

        txBitsE = de2bi(tx_E,2,'left-msb');
        rxBitsE = de2bi(rx_E,2,'left-msb');

        QBER_rad(ii) = ...
            sum(txBitsB(:)~=rxBitsB(:))/numel(txBitsB);

        BER_eve_rad(ii) = ...
            sum(txBitsE(:)~=rxBitsE(:))/numel(txBitsE);

        SecretFraction_rad(ii) = max(0,I_AB-I_AE);

        Pacc_rad(ii) = accepted_rad_B(ii)/Nsym;

        ESKR_rad_kbps(ii) = ...
            (fs/sps)*Pacc_rad(ii)*SecretFraction_rad(ii)/1000;

    end

    %% Rectangular

    valid_rect = ...
        Ierr<r_decision & ...
        Qerr<r_decision;

    accepted_rect_B(ii)=sum(valid_rect);

    if any(valid_rect)

        tx_B = dataBits_B(valid_rect);
        rx_B = rxDecoded_B(valid_rect);

        I_AB = MutualInfo(tx_B,rx_B,M);

        tx_E = dataBits_B(valid_rect);
        rx_E = rxDecoded_E(valid_rect);

        I_AE = MutualInfo(tx_E,rx_E,M);

        ser_rect_B(ii) = mean(tx_B~=rx_B);
        SER_eve_rect(ii) = mean(tx_E~=rx_E);

        txBitsB = de2bi(tx_B,2,'left-msb');
        rxBitsB = de2bi(rx_B,2,'left-msb');

        txBitsE = de2bi(tx_E,2,'left-msb');
        rxBitsE = de2bi(rx_E,2,'left-msb');

        QBER_rect(ii) = ...
            sum(txBitsB(:)~=rxBitsB(:))/numel(txBitsB);

        BER_eve_rect(ii) = ...
            sum(txBitsE(:)~=rxBitsE(:))/numel(txBitsE);

        SecretFraction_rect(ii) = max(0,I_AB-I_AE);

        Pacc_rect(ii) = accepted_rect_B(ii)/Nsym;

        ESKR_rect_kbps(ii) = ...
            (fs/sps)*Pacc_rect(ii)*SecretFraction_rect(ii)/1000;

    end

    %% Angular

    r_min=sigma_B*rho;

    terminate = false;
    
    for it = 1:Nt

        theta = thetaVec(it);
        theta_th = deg2rad(theta);
        valid_ang=...
            abs(phaseErr1)<theta_th & ...
            mag>r_min;
        
        accepted_ang_B(ii,it)=sum(valid_ang);
       
        Nacc = sum(valid_ang);

        if any(valid_ang)

            tx_B = dataBits_B(valid_ang);
            rx_B = rxDecoded_B(valid_ang);

            I_AB = MutualInfo(tx_B,rx_B,M);

            tx_E = dataBits_B(valid_ang);
            rx_E = rxDecoded_E(valid_ang);

            I_AE = MutualInfo(tx_E,rx_E,M);

            ser_ang_B(ii,it) = mean(tx_B~=rx_B);
            SER_eve_ang(ii,it) = mean(tx_E~=rx_E);

            txBitsB = de2bi(tx_B,2,'left-msb');
            rxBitsB = de2bi(rx_B,2,'left-msb');

            txBitsE = de2bi(tx_E,2,'left-msb');
            rxBitsE = de2bi(rx_E,2,'left-msb');

            QBER_ang(ii,it) = ...
                sum(txBitsB(:)~=rxBitsB(:))/numel(txBitsB);

            BER_eve_ang(ii,it) = ...
                sum(txBitsE(:)~=rxBitsE(:))/numel(txBitsE);

            SecretFraction_ang(ii,it) = max(0,I_AB-I_AE);

            Pacc_ang(ii,it) = accepted_ang_B(ii,it)/Nsym;

            ESKR_ang_kbps(ii,it) = ...
                (fs/sps)*Pacc_ang(ii,it)*SecretFraction_ang(ii,it)/1000;

        end
    
    
        if Nacc == 0
            fprintf('\nNo accepted symbols at rho = %.2f\n', rho);
            fprintf('Stopping sweep...\n');

            terminate = true;
        break;
        end
    end

        if terminate
            QBER_ang  = QBER_ang(1:ii-1,:);
            QBER_rad  = QBER_rad(1:ii-1);
            QBER_rect  = QBER_rect(1:ii-1);

            ESKR_rad_kbps = ESKR_rad_kbps(1:ii-1);
            ESKR_rect_kbps = ESKR_rect_kbps(1:ii-1);
            ESKR_ang_kbps = ESKR_ang_kbps(1:ii-1,:);

            rhoVec    = rhoVec(1:ii-1);

            ser_rad_B  = ser_rad_B(1:ii-1);
            ser_rect_B = ser_rect_B(1:ii-1);
            ser_ang_B  = ser_ang_B(1:ii-1,:);

            BER_eve_rad  = BER_eve_rad(1:ii-1);
            BER_eve_rect = BER_eve_rect(1:ii-1);
            BER_eve_ang  = BER_eve_ang(1:ii-1,:);

            SecretFraction_rad = SecretFraction_rad(1:ii-1);
            SecretFraction_rect = SecretFraction_rect(1:ii-1);
            SecretFraction_ang = SecretFraction_ang(1:ii-1,:);

            SER_eve_rad = SER_eve_rad(1:ii-1);
            SER_eve_rect = SER_eve_rect(1:ii-1);
            SER_eve_ang = SER_eve_ang(1:ii-1,:);

            accepted_rad_B = accepted_rad_B(1:ii-1);
            accepted_rect_B = accepted_rect_B(1:ii-1);
            accepted_ang_B = accepted_ang_B(1:ii-1,:);

            Pacc_rad = Pacc_rad(1:ii-1);
            Pacc_rect = Pacc_rect(1:ii-1);
            Pacc_ang = Pacc_ang(1:ii-1,:);
            break;
        end
end
    
%%
figure;

subplot(2,1,1); hold on;

sym0 = (dataBits_B == 0);
sym1 = (dataBits_B == 1);
sym2 = (dataBits_B == 2);
sym3 = (dataBits_B == 3);

plot(real(rxData_B(sym0)),imag(rxData_B(sym0)),'.',...
    'Color',[0 0.4470 0.7410],...
    'MarkerSize',5);

plot(real(rxData_B(sym1)),imag(rxData_B(sym1)),'.',...
    'Color',[0.8500 0.3250 0.0980],...
    'MarkerSize',5);

plot(real(rxData_B(sym2)),imag(rxData_B(sym2)),'.',...
    'Color',[0.4660 0.6740 0.1880],...
    'MarkerSize',5);

plot(real(rxData_B(sym3)),imag(rxData_B(sym3)),'.',...
    'Color',[0.4940 0.1840 0.5560],...
    'MarkerSize',5);

% Ideal constellation
refConst = mod_depth*qammod(0:M-1,M,'UnitAveragePower',true);

plot(real(refConst),imag(refConst),...
    'kx','MarkerSize',12,'LineWidth',2);

axis equal;
grid on;
xlabel('In-Phase');
ylabel('Quadrature');

legend(...
    'Symbol 0',...
    'Symbol 1',...
    'Symbol 2',...
    'Symbol 3',...
    'Ideal',...
    'Location','best');

title('Received Equalized Constellation (Bob)');


subplot(2,1,2); hold on;

sym0 = (dataBits_B == 0);
sym1 = (dataBits_B == 1);
sym2 = (dataBits_B == 2);
sym3 = (dataBits_B == 3);

plot(real(rxData_E(sym0)),imag(rxData_E(sym0)),'.',...
    'Color',[0 0.4470 0.7410],...
    'MarkerSize',5);

plot(real(rxData_E(sym1)),imag(rxData_E(sym1)),'.',...
    'Color',[0.8500 0.3250 0.0980],...
    'MarkerSize',5);

plot(real(rxData_E(sym2)),imag(rxData_E(sym2)),'.',...
    'Color',[0.4660 0.6740 0.1880],...
    'MarkerSize',5);

plot(real(rxData_E(sym3)),imag(rxData_E(sym3)),'.',...
    'Color',[0.4940 0.1840 0.5560],...
    'MarkerSize',5);

% Ideal constellation
refConst = mod_depth*qammod(0:M-1,M,'UnitAveragePower',true);

plot(real(refConst),imag(refConst),...
    'kx','MarkerSize',12,'LineWidth',2);

axis equal;
grid on;
xlabel('In-Phase');
ylabel('Quadrature');

legend(...
    'Symbol 0',...
    'Symbol 1',...
    'Symbol 2',...
    'Symbol 3',...
    'Ideal',...
    'Location','best');

title('Received Equalized Constellation (Eve)');

%% 3-D Plots for optimal rho and theta selection and analysis

thetaVec = 1:1:45;
[RHO, THETA] = meshgrid(rhoVec, thetaVec);

% Psift
figure;

surf(RHO,THETA,Pacc_ang.');

shading interp;
grid on;
colorbar;

xlabel('\rho');
ylabel('\theta (deg)');
zlabel('P_{sift}');

title('Angular Sifting Probability');
view(135,30);


% QBER
figure;

surf(RHO,THETA,QBER_ang.');

shading interp;
grid on;
colorbar;

xlabel('\rho');
ylabel('\theta (deg)');
zlabel('QBER');

title('Angular QBER');
view(135,30);

% ESKR
figure;

surf(RHO,THETA,ESKR_ang_kbps.');

shading interp;
grid on;
colorbar;

xlabel('\rho');
ylabel('\theta (deg)');
zlabel('ESKR (kbps)');

title('Angular Estimated Secret Key Rate');
view(135,30);

% Secret Fraction
figure;

surf(RHO,THETA,SecretFraction_ang.');

shading interp;
grid on;
colorbar;

xlabel('\rho');
ylabel('\theta (deg)');
zlabel('Secret Fraction');

title('Angular Secret Fraction');
view(135,30);