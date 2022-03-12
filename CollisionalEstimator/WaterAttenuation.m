function [WaterMassAttenuationValue, WaterMassAbsorptionValue] = WaterAttenuation(E)
    %Mass Attenuation Value Water from NIST
    load('WaterMassAtten.mat')
    q = length(WaterMassAtten(:,1));
    IndexRegister = 0;
    for k = 1 : q
        if E > WaterMassAtten(k,1)
        IndexRegister = k;
        end
    end
    BottomMu = WaterMassAtten(IndexRegister,2);
    TopMu = WaterMassAtten(IndexRegister + 1, 2);

    BottomMuAbs = WaterMassAtten(IndexRegister, 3);
    TopMuAbs = WaterMassAtten(IndexRegister + 1, 3);

    BottomEnergy = WaterMassAtten(IndexRegister, 1);
    TopEnergy = WaterMassAtten(IndexRegister + 1, 1);

    MuSlope = (TopMu - BottomMu)/(TopEnergy - BottomEnergy);
    WaterMassAttenuationValue = MuSlope * (E - BottomEnergy) + BottomMu;

    MuSlopeAbs = (TopMuAbs - BottomMuAbs)/(TopEnergy - BottomEnergy);
    WaterMassAbsorptionValue = MuSlopeAbs * (E - BottomEnergy) + BottomMuAbs;
end
