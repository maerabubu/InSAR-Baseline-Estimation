function XYZ = lp2xyz(p_orbit,a_time, r_time,approxXYZCentre)

if nargin < 5
    h = zeros(size(line));
end

if nargin < 6
    ell.a = 6378137.0 + h;			% semimajor axis wgs84
    ell.b = 6356752.3142451794975639665996337 + h;			% semiminor axis wgs84
end

if nargin < 7
    CRITERPOS = 1e-6;
end

if nargin < 8
    MAXITER = 100;
end
[row, col] = size(a_time);
equationset = zeros(3, 1, row, col);
solxyz      = zeros(3, row, col);


VLIGHT = 299792458;             
pos = getXYZ( a_time,p_orbit); 
%px=a_time.^3*Coeff.X(1)+a_time.^2*Coeff.X(2)+a_time*Coeff.X(3)+Coeff.X(4);
%py=a_time.^3*Coeff.Y(1)+a_time.^2*Coeff.Y(2)+a_time*Coeff.Y(3)+Coeff.Y(4);
%pz=a_time.^3*Coeff.Z(1)+a_time.^2*Coeff.Z(2)+a_time*Coeff.Z(3)+Coeff.Z(4);
posv = getXYZDot(a_time,p_orbit);
%vx=3*Coeff.X(1).*a_time.^2+2*Coeff.X(2).*a_time+Coeff.X(3);
%vy=3*Coeff.Y(1).*a_time.^2+2*Coeff.Y(2).*a_time+Coeff.Y(3);
%vz=3*Coeff.Z(1).*a_time.^2+2*Coeff.Z(2).*a_time+Coeff.Z(3);
    %x,y,z of image center    
    X = ones(row, col) * approxXYZCentre.X;
    Y = ones(row, col) * approxXYZCentre.Y;
    Z = ones(row, col) * approxXYZCentre.Z;

    
    partialsxyz = zeros(3, 3, row, col);
	
    for iter = 1:MAXITER
    	dsat_P.x= X - pos.X;
        dsat_P.y= Y - pos.Y;
        dsat_P.z= Z - pos.Z;
        
        
		equationset(1,1,:,:) = (posv.X.*dsat_P.x + posv.Y.*dsat_P.y+ posv.Z.*dsat_P.z).*(-1);
		equationset(2,1,:,:) = (dsat_P.x.*dsat_P.x + dsat_P.y.*dsat_P.y + dsat_P.z.*dsat_P.z - (VLIGHT .*r_time).*(VLIGHT .*r_time)).*(-1);
		equationset(3,1,:,:) = ((X.*X + Y.*Y)./(ell.a.*ell.a) + (Z./ell.b).*(Z./ell.b) -1)*(-1);
        
		
        partialsxyz(1,1,:,:) = posv.X;
		partialsxyz(1,2,:,:) = posv.Y;
		partialsxyz(1,3,:,:) = posv.Z;
		partialsxyz(2,1,:,:) = 2.*dsat_P.x;
		partialsxyz(2,2,:,:) = 2.*dsat_P.y;
		partialsxyz(2,3,:,:) = 2.*dsat_P.z;
		partialsxyz(3,1,:,:) = (2*X)./(ell.a.*ell.a);
		partialsxyz(3,2,:,:) = (2*Y)./(ell.a.*ell.a);
		partialsxyz(3,3,:,:) = (2*Z)./(ell.b.*ell.b);

		
        %get inv(partialsxyz) for loop maybe OK? try first
        for i=1:row
            for j=1:col
                solxyz(:,i,j) = inv(partialsxyz(:,:,i,j)) * equationset(:,:,i,j);
            end
        end
 	
		X = shiftdim(solxyz(1,:,:)) + X;	
		Y = shiftdim(solxyz(2,:,:)) + Y;	
		Z = shiftdim(solxyz(3,:,:)) + Z;	
        
	
        if abs(solxyz) < CRITERPOS
            break;
        end
        if iter > 10
        %disp(line, pixel);
        disp('reached the max iteration number in lp2xyz()');
        break;
        end
    end
    XYZ.X=X;XYZ.Y=Y;
    XYZ.Z=Z;
    end %for iter = 1 :MAXITER
    %check iter numbe
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
    

