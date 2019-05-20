function [B, Bpar, Bperp, theta] = BBparBperpTheta(Master, Point, Slave)
B = norm(Master - Slave); %Length of Baseline

range1 = norm(Master - Point);
range2 = norm(Slave - Point);

Bpar = range1 - range2;  %// parallel baseline, sign ok

r1 = Master - Point;
r2 = Slave - Point;

theta = vec_angle(Point, r1);

Bperp = B^2 - Bpar^2;

if Bperp < 0.0
    Bperp = 0.;
elseif theta > vec_angle(Point, r2)
    Bperp = sqrt(Bperp);
else
    Bperp = 0 - sqrt(Bperp);
end
end