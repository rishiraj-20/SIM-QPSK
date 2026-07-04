%% CHECK PHASE DRIFT ACROSS DATA REGION

% Split data into 10 chunks, measure phase of each
Nchunk = 10;
chunkSize = floor(Nsym / Nchunk);

fprintf('--- Phase drift across data region ---\n');
for k = 1:Nchunk
    chunk = rxData((k-1)*chunkSize+1 : k*chunkSize);
    distC = abs(chunk - refConst.');
    [~, idxC] = min(distC, [], 2);
    phi_k = angle(mean(chunk .* conj(refConst(idxC))));
    fprintf('Chunk %2d : phase = %+.2f deg\n', k, phi_k*180/pi);
end

% FIX 1 — Skip first chunk (transient region)       
% In most of the cases, first chunk was more noisy

skipSym = chunkSize;   % discarding first ~10% symbols
rxData     = rxData(skipSym+1:end);
dataBits   = dataBits(skipSym+1:end);
Nsym       = length(rxData);

% FIX 2 — Correct the steady-state offset
% Estimate from chunks 2-10 only (stable region)
stableChunks = rxData;   % already skipped chunk 1
distS = abs(stableChunks - refConst.');
[~, idxS] = min(distS, [], 2);
phi_steady = angle(mean(stableChunks .* conj(refConst(idxS))));
fprintf('Steady state phase offset = %.2f deg\n', phi_steady*180/pi);

rxData = rxData .* exp(-1j * phi_steady);