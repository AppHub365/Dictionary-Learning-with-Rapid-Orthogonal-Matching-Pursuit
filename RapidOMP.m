%% Rapid Orthogonal Matching Pursuit Estimate
% Authored by Ayan Chatterjee (ayan@outlook.com) on 23 December 2019

function sparseRep = RapidOMP(signals, basis, maxBasisLim, maxErrorTol)

%% pre-process
[~, numSignals] = size(signals);
[~, numBasis] = size(basis);
basisMatch = single(zeros(numBasis, numSignals));
sparseRep = single(zeros(numBasis, numSignals));
residue = signals;
if(~exist('maxErrorTol', 'var'))
    maxErrorTol = 1e-05; % change as required
end
signalsRemaining = 1:numSignals;

%% starting rapid OMP algorithm
for lim = 1:maxBasisLim
    %% basis selection
    basisProject = abs(basis' * residue);
    for basisCount = 1:numBasis
        selectedBasisPos = sum(basisProject > basisProject(basisCount, :)) == 0;
        basisMatch(basisCount, signalsRemaining(selectedBasisPos == 1)) = 1; % these samples will get allocated #basisCount basis
    end
    
    %% estimate residue
    [uniqueBasisMatch, ~, invBasisPos] = unique(basisMatch', 'rows');
    for invIter = 1:size(uniqueBasisMatch, 1)
        matchPos = find(uniqueBasisMatch(invIter, :) == 1);
        invBasis = basis(:, matchPos);
        invBasis = pinv(invBasis' * invBasis) * invBasis';
        sparseRep(matchPos, invBasisPos == invIter) = invBasis * signals(:, invBasisPos == invIter);
    end
    residue = signals - basis * sparseRep; % residue (r) = y - Da
    signalsRemaining = find(sum(abs(residue)) > maxErrorTol);
    residue = residue(:, signalsRemaining);
end

end