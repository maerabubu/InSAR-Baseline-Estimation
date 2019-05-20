function [time,solution] = xyz2t(pointOnEllips,info,coeff)
if nargin < 7
    CRITERTIM = 10e-10;
end
SOL=299792458;
if nargin < 6
    MAXITER = 20;
end
[m,n]=size(pointOnEllips.X);
timeAzimuth =ones(m,n).*line2ta(info.sensor_para.current_window.bottom/2,info.sensor_para.t_Azi1,info.sensor_para.line_time_interval);      
for iter = 1:MAXITER
            satellitePosition = getXYZ(timeAzimuth,coeff);
            satelliteVelocity = getXYZDot(timeAzimuth,coeff);
            satelliteAcceleration = getXYZDotDot(timeAzimuth,coeff);
            delta.X = pointOnEllips.X-satellitePosition.X;
            delta.Y = pointOnEllips.Y-satellitePosition.Y;
            delta.Z = pointOnEllips.Z-satellitePosition.Z;
            solution{iter} = -eq1_Doppler(satelliteVelocity, delta) ./ eq1_Doppler_dt(delta, satelliteVelocity, satelliteAcceleration);
            timeAzimuth =timeAzimuth+ solution{iter};
            if max(max(abs(solution{iter}))) < CRITERTIM
                break;
            end
            if (iter >= MAXITER) 
            disp('ERROR');
                break;
            end
end
        
satellitePosition = getXYZ(timeAzimuth,coeff);
delta.X = pointOnEllips.X-satellitePosition.X;
delta.Y = pointOnEllips.Y-satellitePosition.Y;
delta.Z = pointOnEllips.Z-satellitePosition.Z;
timeRange = sqrt(delta.X.^2 +delta.Y.^2+ delta.Z.^2)./ SOL;
time.timeAzimuth=timeAzimuth;
time.timeRange=timeRange;
end
function [POS]=getXYZ(azTime,coeff)
azTimeNormal = (azTime - coeff.time(floor(length(coeff.time) / 2)+1))./10;
%azTimeNormal = azTime;
POS.X=polyVal1D(azTimeNormal, coeff.X);
POS.Y=polyVal1D(azTimeNormal, coeff.Y);
POS.Z=polyVal1D(azTimeNormal, coeff.Z);
end
function [sum]=polyVal1D(aziTime, coeffs)
sum = 0.0;
d=length(coeffs);
        while d > 0
            sum =sum.* aziTime;
            sum =sum+ coeffs(d);
            d=d-1;
        end
end
function [POS]=getXYZDot(azTime,coeff)
azTimeNormal = (azTime - coeff.time(floor(length(coeff.time) / 2)+1))./10 ;
%azTimeNormal = azTime;
DEGREE = length(coeff.X)-1;
x = coeff.X(2);
y = coeff.Y(2);
z = coeff.Z(2);
        for  i = 2:DEGREE
            powT = i .*power(azTimeNormal, i - 1);
            x = x+coeff.X(i+1) .* powT;
            y = y+coeff.Y(i+1) .* powT;
            z = z+coeff.Z(i+1) .* powT;
        end
POS.X=x/10;
POS.Y=y/10; 
POS.Z=z/10;
end
function  [POS]=getXYZDotDot(azTime,coeff)
poly_degree = length(coeff.X)-1;
azTimeNormal = (azTime - coeff.time(floor(length(coeff.time) / 2))+1)./10 ;
%azTimeNormal = azTime;
x=0; y=0; z=0;
        for  i = 2: poly_degree
            powT = ((i - 1) * i) .* power(azTimeNormal, i - 2);
            x =x+ coeff.X(i+1) .* powT;
            y =y+ coeff.Y(i+1) .* powT;
            z =z+ coeff.Z(i+1) .* powT;
        end
POS.X=x/100;
POS.Y=y/100;
POS.Z=z/100;
end
function [POS]=eq1_Doppler(satVelocity,pointOnEllips)
POS=satVelocity.X.*pointOnEllips.X+satVelocity.Y.*pointOnEllips.Y+satVelocity.Z.*pointOnEllips.Z;
end
function [POS]=eq1_Doppler_dt(pointEllipsSat,satVelocity,satAcceleration)
     POS=  satAcceleration.X.*pointEllipsSat.X+satAcceleration.Y.*pointEllipsSat.Y+satAcceleration.Z.*pointEllipsSat.Z ...
       - satVelocity.X.*satVelocity.X - satVelocity.Y.*satVelocity.Y - satVelocity.Z.*satVelocity.Z;
end