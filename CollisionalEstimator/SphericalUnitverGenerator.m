function [x1, y1, z1, Theta] = SphericalUnitverGenerator
Theta = RandomPi;
Phi = rand * 2 * pi;
x1 = cos(Theta);
y1 = sin(Theta) * cos(Phi);
z1 = sin(Theta) * sin(Phi);
end
