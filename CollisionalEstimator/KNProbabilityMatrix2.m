function [KNDivMatrix, DividerIndex, thetaArray, KNCummulativeMatrix, KNArrayMatrix, EnergyMatrix] = KNProbabilityMatrix2(n,E,dividerVal, delta)
NSequence = E / delta;
EnergyMatrix = zeros(NSequence, 1);
KNDivMatrix = zeros(NSequence, dividerVal);
KNCummulativeMatrix = zeros(NSequence, n);
KNArrayMatrix = zeros(NSequence, n);
thetaArrayMatrix = zeros(NSequence, n);

    for Sequence = 1 : NSequence
        EnergyMatrix(Sequence,:) = E;
        [DividerIndex, KNCummulative, DividerKN, thetaArray, KNArray] = KNProbability(n,E,dividerVal);
        KNDivMatrix(Sequence,:) = DividerKN;
        thetaArrayMatrix(Sequence,:) = thetaArray;
        KNCummulativeMatrix(Sequence,:) = KNCummulative;
        KNArrayMatrix(Sequence,:) = KNArray;
        E = E - delta;
    end

end