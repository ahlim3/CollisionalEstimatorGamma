function CollisionalEstimatorWater(radius, E, N)
n = 10000;
dividerVal = 10;
delta = 0.01;
NEnergyStep = E/delta;

MaximumLoopNumber = 10;
xcoordinateParticle = NaN(N, MaximumLoopNumber);
ycoordinateParticle = NaN(N, MaximumLoopNumber);
zcoordinateParticle = NaN(N, MaximumLoopNumber);
dxMatrix = NaN(N, MaximumLoopNumber);
dyMatrix = NaN(N, MaximumLoopNumber);
dzMatrix = NaN(N, MaximumLoopNumber);
thetaMatrix = NaN(N, MaximumLoopNumber);
phiMatrix = NaN(N, MaximumLoopNumber);
LengthMatrix = NaN(N, MaximumLoopNumber);
d1Matrix = NaN(N, MaximumLoopNumber);
ExceedPath = NaN(N, MaximumLoopNumber);
PathLengthMatrix = zeros(N, MaximumLoopNumber);
PathTrack = zeros(N, 1);
EndingEnergy = NaN(N, 1);
EndingCoordinateX = NaN(N, MaximumLoopNumber);
EndingCoordinateY = NaN(N, MaximumLoopNumber);
EndingCoordinateZ = NaN(N, MaximumLoopNumber);
EnergyInteraction = zeros(N, MaximumLoopNumber);
EnergyTrack = zeros(N, 1);
InteractionPerPart = zeros(N, 1);
DFromOrigin = NaN(N, MaximumLoopNumber);
FluenceDivider = (radius / 0.25);
FluenceStep = radius / FluenceDivider;
FluenceMatrix = zeros(1,FluenceDivider + 1);
FluenceMatrix2 = zeros(1,FluenceDivider + 1);
FleunceRadius = 0:FluenceStep:radius;
FleunceArea = zeros(1,length(FleunceRadius));
ELoss = zeros(N, MaximumLoopNumber);

for i = 1 : length(FleunceRadius)
    FleunceArea(i) = 4 * pi * FleunceRadius(i)^2 * FluenceStep;
end

FluenceIndex = zeros(1,20);

for i = 1 : 20
    FluenceIndex(i) = i * FluenceStep * FluenceDivider / 20;
end

[KNDivMatrix, DividerIndex, thetaArrayMatrix, KNCummulativeMatrix, ~, EnergyMatrix] = KNProbabilityMatrix2(n,E,dividerVal, delta);




for k = 1 : N
    %Initial Condition
    irr = 1;
    IntermediateE = E;
    %Setting the starting point of particle
    [x, y, z] = SphereInsidePoint(radius);
    xcoordinateParticle(k,irr) = x;
    ycoordinateParticle(k,irr) = y;
    zcoordinateParticle(k,irr) = z;
    %Registering the first point
    EndingCoordinateX(k,irr) = xcoordinateParticle(k,irr);
    EndingCoordinateY(k,irr) = ycoordinateParticle(k,irr);
    EndingCoordinateZ(k,irr) = zcoordinateParticle(k,irr);
    phiIntermediate = rand * 2 * pi;
    thetaIntermediate = RandomPi;

    %Unit Vector for up to 1st collision site
    dx = cos(thetaIntermediate);
    dy = sin(thetaIntermediate) * cos(phiIntermediate);
    dz = sin(thetaIntermediate) * sin(phiIntermediate);
    dxMatrix(k,irr) = dx;
    dyMatrix(k,irr) = dy;
    dzMatrix(k,irr) = dz;
    
    %Condition from the origin
    distFromOrigin = sqrt(xcoordinateParticle(k,irr)^2 + ycoordinateParticle(k,irr)^2 + zcoordinateParticle(k,irr)^2);
    LengthMatrix(k,irr) = distFromOrigin;
    thetaMatrix(k, irr) = thetaIntermediate;
    phiMatrix(k,irr) = phiIntermediate;
    EnergyInteraction(k,irr) = IntermediateE;

    LoopNumber = 2;
    %Condition for particle termination
    while IntermediateE > delta
        if LoopNumber > MaximumLoopNumber
            break
        end
        %Collision Site # N + 1;
        irr = irr + 1;
        %Absoprtion and Attenuation Coefficient as 1.0 g/cm^3 of water
        [mu, ~] = WaterAttenuation(IntermediateE);
        %Average distance from the original site to next collisional site
        d1 = 1/mu;
        %Actual collisional distance
        d1Matrix(k,irr) = d1 * log (1/rand);
        
        EnergyReturn = 1;

        for l = 1 : NEnergyStep
            if IntermediateE > EnergyMatrix(l)
                EnergyReturn = l;
            end
        end


        %theta selection for Compton (KN)
        compare = rand;

        for l = 1 : dividerVal - 1
            if compare > KNDivMatrix(EnergyReturn,l)
                startPoint = DividerIndex(l);
            end
        end
        for m = startPoint:n
            if compare < KNCummulativeMatrix(EnergyReturn,m)
                thetaReturn = m;
                break
            end
        end
        
        thetaMatrix(k, irr) = thetaArrayMatrix(thetaReturn);
        thetaPrevious = thetaMatrix(k, irr - 1);
        thetaIntermediate = thetaMatrix(k, irr);
        phiIntermediate = rand * 2 * pi;
        phiPrevious = phiMatrix(k, irr - 1);
        phiMatrix(k,irr) = phiIntermediate;
        
        %Reference frame in terms of unit vector from previous condition
        u0 = cos(thetaPrevious);
        v0 = sin(thetaPrevious) * cos(phiPrevious);
        w0 = sin(thetaPrevious) * sin(phiPrevious);
        
        CosIntermediateTheta = cos(thetaIntermediate);
        SinIntermediateTheta = sin(thetaIntermediate);
        CosIntermdiatePhi = cos(phiIntermediate);
        SinIntermediatePhi = sin(phiIntermediate);
        w0sqrt = sqrt(1-w0*w0);
        
        %Conversion of reference frame to x,y,z frame
        dx = u0 * CosIntermediateTheta + (u0 * w0 * CosIntermdiatePhi - v0 * SinIntermediateTheta) * (SinIntermediatePhi / w0sqrt);
        dy = v0 * CosIntermediateTheta + (v0 * w0 * CosIntermdiatePhi + u0 * SinIntermediatePhi) * (SinIntermediatePhi / w0sqrt);
        dz = w0 * CosIntermediateTheta - CosIntermdiatePhi * SinIntermediateTheta * w0sqrt;
        
        %Registering condition in terms of x,y,z frame
        dxMatrix(k,irr) = dx;
        dyMatrix(k,irr) = dy;
        dzMatrix(k,irr) = dz;
        
        %Particle Location at the Collisional Site
        xcoordinateParticle(k,irr) = xcoordinateParticle(k,irr - 1) + d1 * dx;
        ycoordinateParticle(k,irr) = ycoordinateParticle(k,irr - 1) + d1 * dy;
        zcoordinateParticle(k,irr) = zcoordinateParticle(k,irr - 1) + d1 * dz;
        EndingCoordinateX(k,irr) = xcoordinateParticle(k,irr);
        EndingCoordinateY(k,irr) = ycoordinateParticle(k,irr);
        EndingCoordinateZ(k,irr) = zcoordinateParticle(k,irr);
        PathLengthMatrix(k,irr - 1) = sqrt((xcoordinateParticle(k,irr) - xcoordinateParticle(k,irr - 1))^2 + (ycoordinateParticle(k,irr) - ycoordinateParticle(k,irr - 1))^2 + (zcoordinateParticle(k,irr) - zcoordinateParticle(k,irr - 1))^2);

        IntermediateE = ComptonScatteringHv(IntermediateE,thetaIntermediate);
        
        distFromOrigin = sqrt(xcoordinateParticle(k,irr)^2 + ycoordinateParticle(k,irr)^2 + zcoordinateParticle(k,irr)^2);
        LengthMatrix(k,irr) = distFromOrigin;
       
        if distFromOrigin > radius
            dxCompare = dxMatrix(k,irr);
            dyCompare = dyMatrix(k, irr);
            dzCompare = dzMatrix(k, irr);
            xCompare = xcoordinateParticle(k,irr - 1);
            yCompare = ycoordinateParticle(k, irr - 1);
            zCompare = zcoordinateParticle(k, irr - 1);
            %Distance from origin to P1
            P1magnitude = sqrt(xCompare^2 + yCompare^2 + zCompare^2);
            %Unitvector mutation correction
            dMagnitude = sqrt(dxCompare^2 + dyCompare^2 + dzCompare^2);

            %Origin to P1 interms of Unitvector
            xUnitvector = xCompare / P1magnitude;
            yUnitvector = yCompare / P1magnitude;
            zUnitVector = zCompare / P1magnitude;
            dxCompare = dxCompare / dMagnitude;
            dyCompare = dyCompare / dMagnitude;
            dzCompare = dzCompare / dMagnitude;

            %Angle between two vectors, origin to P1 and P1 to intersection
            %Realizaing angle is paralleogram pi - angle between vectors
            IntermediateAngleP1 = pi - acos(xUnitvector * dxCompare + yUnitvector * dyCompare + zUnitVector * dzCompare);
            
            %d is more than radius use of SSA
            SineRatioRadius = sin(IntermediateAngleP1) / radius;
            AngleP1ToSurface = asin(SineRatioRadius * P1magnitude);
            AngleOriginToP1 = pi - IntermediateAngleP1 - AngleP1ToSurface;
            dFromP1ToSurface = sin(AngleOriginToP1) / SineRatioRadius;

            ExceedPath(k,1) = dFromP1ToSurface;
            PathLengthMatrix(k,irr - 1) = dFromP1ToSurface;
            xcoordinateParticle(k,irr) = xCompare + dxCompare * dFromP1ToSurface;
            ycoordinateParticle(k,irr) = yCompare + dyCompare * dFromP1ToSurface;
            zcoordinateParticle(k,irr) = zCompare + dzCompare * dFromP1ToSurface;
            EnergyInteraction(k,irr) = IntermediateE;
            ELoss(k,irr - 1) = EnergyInteraction(k, irr-1) - IntermediateE;
            break
        end
        DFromOrigin(k,irr-1) = sqrt(xcoordinateParticle(k,irr)^2 + ycoordinateParticle(k,irr)^2 + zcoordinateParticle(k,irr)^2);
        DChecking = DFromOrigin(k,irr-1);
        for i = 1 : FluenceDivider
            if DChecking > FleunceRadius(i)
                QLDED = i;
            end
        end
        FluenceMatrix(1, QLDED) = FluenceMatrix(1, QLDED) + 1;
        LoopNumber = LoopNumber + 1;
        %Energy of particle after the collision
        EnergyInteraction(k,irr) = IntermediateE;
        ELoss(k,irr - 1) = EnergyInteraction(k, irr-1) - IntermediateE;
    end
        EndingEnergy (k) = IntermediateE;
        EnergyTrack(k) = sum(ELoss(k,:));
        PathTrack(k) = sum(PathLengthMatrix(k,:));
        if DFromOrigin(k, irr-1) > radius
            DFromOrigin(k, irr-1) = NaN;
        end
        InteractionPerPart(k) = irr - 1;
        
end

xr = zeros(1000,1);
yr = zeros(1000,1);
zr = zeros(1000,1);
for l = 1:1000
    [x, y, z] = SphereGenerator(radius);
    xr(l) = x;
    yr(l) = y;
    zr(l) = z;
end

FluenceMatrix2(1) = FluenceMatrix(1) / (FleunceArea(i));
for i = 2 : length(FleunceRadius)
    FluenceMatrix2(i) = FluenceMatrix(i) / (FleunceArea(i) - FleunceArea(i-1));
end
FluenceMatrix2(1) = [];
FleunceRadius(1) = [];
FluenceMatrix2(end) = [];
FleunceRadius(end) = [];

FluenceMatrix2 = FluenceMatrix2';
FleunceRadius = FleunceRadius';
xcoordinateParticleCollision = xcoordinateParticle;
ycoordinateParticleCollision = ycoordinateParticle;
zcoordinateParticleCollision = zcoordinateParticle;
xcoordinateParticleCollision(:,1) = [];
ycoordinateParticleCollision(:,1) = [];
zcoordinateParticleCollision(:,1) = [];


ax = gca; 
ax.FontSize = 16;

figure(1)
plot3(xcoordinateParticleCollision,ycoordinateParticleCollision,zcoordinateParticleCollision,'.',xr,yr,zr,'-','MarkerSize',15)
legend('Collision 1','Collision 2','Collision 3','Collision 4','Collision 5','Collision 6','Collision 7','Collision 8','Collision 9','Sphere')
PlotTitle = sprintf('Particle Projetion with %.0f Bq with position termination at the Surface', N);
title(PlotTitle)
xlabel('x (cm)')
ylabel('y (cm)')
zlabel('z (cm)')
ax = gca; 
ax.FontSize = 16;



figure(2)
EnergyMu = mean(EnergyTrack);
EnergySigma = std(EnergyTrack);
histogram(EnergyTrack)
xlabel('Energy Deposition Per Particle (MeV)', 'FontSize', 15);
ylabel('Frequency', 'FontSize', 15);
xline(EnergyMu, 'color', 'r', "LineWidth", 2);
xline(EnergyMu - EnergySigma, 'color', 'g', 'LineWidth', 2, 'LineStyle', '--');
xline(EnergyMu + EnergySigma, 'color', 'b', 'LineWidth', 2, 'LineStyle', '--');
ETrackMean= sprintf('Mean = %.3f\n', EnergyMu);
ETrackUpperBound= sprintf('Lower Bound = %.3f\n', EnergyMu - EnergySigma);
ETrackLowerBound= sprintf('Upper Bound = %.3f\n', EnergyMu + EnergySigma);

ETrackTitle = sprintf('Average energy deposition MeV/Particle with %.0f Bq', N);
title(ETrackTitle)
legend('Data',ETrackMean, ETrackLowerBound, ETrackUpperBound)
ax = gca; 
ax.FontSize = 16;

figure(3)
[p,~,mu] = polyfit(FleunceRadius, FluenceMatrix2,9);
f = polyval(p,FleunceRadius,[],mu);
hold on
plot(FleunceRadius, FluenceMatrix2, 'o')
plot(FleunceRadius, f)
legend('Data', '7th Polynomial Fit')
FluenceTitle = sprintf('Fluence Distribution with %.0f Bq', N);
title(FluenceTitle)
xlabel('Radius (cm)')
ylabel('Fluence in photons/cm^2/s')
ax = gca; 
ax.FontSize = 16;
hold off

figure(4)
PathMu = mean(PathTrack);
PathSigma = std(PathTrack);
histogram(PathTrack)
xlabel('Pathlength Per particle (cm)', 'FontSize', 15);
ylabel('Frequency', 'FontSize', 15);
xline(PathMu, 'color', 'r', "LineWidth", 2);
xline(PathMu - PathSigma, 'color', 'g', 'LineWidth', 2, 'LineStyle', '--');
xline(PathMu + PathSigma, 'color', 'b', 'LineWidth', 2, 'LineStyle', '--');
PathMean= sprintf('Mean = %.3f\n', PathMu);
PathUpperBound= sprintf('Lower Bound = %.3f\n', PathMu - PathSigma);
PathLowerBound= sprintf('Upper Bound = %.3f\n', PathMu + PathSigma);

ETrackTitle = sprintf('Average pathlength (cm) /Particle with %.0f Bq', N);
title(ETrackTitle)
legend('Data',PathMean, PathUpperBound, PathLowerBound)
ax = gca; 
ax.FontSize = 16;

end