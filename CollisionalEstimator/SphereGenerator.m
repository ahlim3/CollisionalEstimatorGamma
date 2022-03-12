function [x1, y1, z1] = SphereGenerator(radius)
Theta = RandomPi;
Phi = rand * 2 * pi;
x1 = radius * cos(Theta);
y1 = radius * sin(Theta) * cos(Phi);
z1 = radius * sin(Theta) * sin(Phi);
end

