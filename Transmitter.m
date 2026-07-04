%% PARAMETERS

mod_depth = 0.1;

fc = 950e6;
fs = 1e6;
sps = 12;
M = 4;

Ndata = 5e4;

rolloff = 0.35;
span = 8;

% Local ZC Sequence

Nsync = 1021;
u = 41;

n = (0:Nsync-1).';

zcSync = exp(-1j*pi*u*n.*(n+1)/Nsync);

% Pilot Data

Npilot = 7500;
rng(123);

pilotBits = randi([0 1], Npilot, 1);
pilotSym  = pskmod(pilotBits,2,0,'gray');

% REGENERATE DATA

rng(12345);

dataBits = randi([0 M-1], Ndata, 1);
dataSym  = mod_depth * qammod(dataBits, M, 'UnitAveragePower', true);


frameLengthSym = (Nsync + Npilot + Ndata);

frameSamples = ...
    (frameLengthSym-1)*sps + span*sps + 1;

% %% Frame

txSym = [zcSync; pilotSym; dataSym];

rrcFilter = rcosdesign(rolloff,span,sps,'sqrt');

txSignal = upfirdn(txSym,rrcFilter,sps,1);

% txSignal = txSignal ./ max(abs(txSignal));

tx = sdrtx('Pluto');

tx.CenterFrequency = fc;

tx.BasebandSampleRate = fs;

tx.Gain = -20;
transmitRepeat(tx,txSignal);