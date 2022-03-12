function DsigmaDomega = KleinNishinaCompton(E, theta)
FEgammatheta = 1/(1+E*(1-cos(theta)));
DsigmaDomega = 0.5 * FEgammatheta^2 * (FEgammatheta + FEgammatheta^-1 - sin(theta)^2);
end
