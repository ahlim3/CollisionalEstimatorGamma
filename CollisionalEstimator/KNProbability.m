function [DividerIndex, KNCummulative, DividerKN, thetaArray, KNArray] = KNProbability(n,E,dividerVal)
    step = pi/n;
    KNCummulative = zeros(1,n);
    thetaArray = zeros(1,n);
    dividerStep = n/dividerVal;
    DividerIndex = zeros(1, dividerVal);
    DividerKN = zeros(1, dividerVal);
    KNArray = zeros(1,n);
        KNatTheta = KleinNishinaCompton(E,step);
        KNArray(1) = KNatTheta;
        KNCummulative(1) = KNCummulative(1);
        thetaArray(1) = step;

    for q = 2 : n
        %dOsigma/dOmega at theta in small step
        KNatTheta = KleinNishinaCompton(E,step*q);
        KNArray(q) = KNatTheta;
        %Integrating via summation, adding dSigma/dOmega * delta from
        %previous integral
        KNCummulative(q) = KNCummulative(q - 1) + KNatTheta*step;
        %Theta
        thetaArray(q) = q*step;
    end

    %Normalize to 1.0
    KNCummulative = KNCummulative/KNCummulative(end);

    for k = 1 : dividerVal
        DividerKN(k) = KNCummulative(dividerStep * k);
        DividerIndex(k) = k*dividerStep;
    end
end