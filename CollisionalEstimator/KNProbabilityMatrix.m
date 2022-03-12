function [KNDivMatrix, DividerIndex, thetaArrayMatrix, KNCummulativeMatrix, KNArrayMatrix, EnergyMatrix] = KNProbabilityMatrix(n,E,dividerVal, delta, ELossFraction)
Sequence = 1;
    while E > delta
        EnergyMatrix(Sequence,:) = E;
        [DividerIndex, KNCummulative, DividerKN, thetaArray, KNArray] = KNProbability(n,E,dividerVal);
        KNDivMatrix(Sequence,:) = DividerKN;
        DivIndMatrix(Sequence,:) = DividerIndex;
        thetaArrayMatrix(Sequence,:) = thetaArray;
        KNCummulativeMatrix(Sequence,:) = KNCummulative;
        KNArrayMatrix(Sequence,:) = KNArray;
        EnergyMatrix(Sequence,:) = E;
        E = (1-ELossFraction)*E;
        Sequence = Sequence + 1;
    end

end