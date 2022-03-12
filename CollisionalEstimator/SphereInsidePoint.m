function [x1, y1, z1] = SphereInsidePoint(radius)
Theta = RandomPi;
Phi = rand * 2 * pi;
radius = rand * radius;
x1 = radius * cos(Theta);
y1 = radius * sin(Theta) * cos(Phi);
z1 = radius * sin(Theta) * sin(Phi);
end

