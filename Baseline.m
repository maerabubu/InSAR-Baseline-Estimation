function p_baseline = Baseline(m_info, s_info)
%p_baseline = baseline(m_info, s_info)
%constructer of Class baseline
% /****************************************************************
% // --- B(l,p,h) = a000 + 
% //                a100*l   + a010*p   + a001*h   +
% //                a110*l*p + a101*l*h + a011*p*h +
% //                a200*l^2 + a020*p^2 + a002*h^2
% ****************************************************************/

ell.a = 6378137.0;			% semimajor axis wgs84
ell.b = 6356752.3142451794975639665996337;			% semiminor axis wgs84

N_coeffs = 10; %==10 degree of 3D-poly to model B(l,p,h)

m_orbit = orbit(m_info,3);
s_orbit = orbit(s_info,3);

%sensor wave length
p_baseline.master_wavelength = m_info.sensor_para.Wavelength;

%current window
current_window = m_info.sensor_para.current_window;
%current_window.top=(burstIndex-1)*m_info.subSwath.linesPerBurst+1;
%current_window.bottom=burstIndex*m_info.subSwath.linesPerBurst;
%current_window.left=1;
%current_window.right=m_info.subSwath.samplesPerBurst;
%// ______ Set min/max for normalization ______
p_baseline.L_min = current_window.top;% also used during polyval
p_baseline.L_max = current_window.bottom;% also used during polyval
p_baseline.P_min = current_window.left;% also used during polyval
p_baseline.P_max = current_window.right;% also used during polyval
p_baseline.H_min = 0.0;% also used during polyval
p_baseline.H_max = 5000.0;% also used during polyval

%// ______ Model r=nearrange+drange_dp*p, p starts at 1 ______
p_baseline.nearrange = pix2range(1.0, m_info.sensor_para.t_Range1, m_info.sensor_para.RSR2x);
p_baseline.drange_dp  = pix2range(2.0, m_info.sensor_para.t_Range1, m_info.sensor_para.RSR2x) ...
    - pix2range(1.0, m_info.sensor_para.t_Range1, m_info.sensor_para.RSR2x);

p_baseline.nearrange = p_baseline.nearrange - p_baseline.drange_dp;% (p starts at 1)
%%%may be something wrong. Please check it when debuging.... 


%// ______ Loop counters ______
cnt = 1;%// matrix index
N_pointsL     = 5; %// every 10km in azimuth
N_pointsP     = 10; %// every 10km ground range
N_heights     = 4; %// one more level than required for poly
deltapixels = (current_window.right -  current_window.left + 1) / N_pointsP;
deltalines  = (current_window.bottom -  current_window.top + 1)  / N_pointsL;
deltaheight = (p_baseline.H_max-p_baseline.H_min) / N_heights;

%// ______ Matrices for modeling Bperp (BK 21-mar-01) ______
%// --- For stability of normalmatrix, fill AMATRIX with normalized line, etc.
BPERP     = zeros(N_pointsL*N_pointsP*N_heights,1); %// perpendicular baseline
BPAR      = zeros(N_pointsL*N_pointsP*N_heights,1); %// parallel baseline
THETA     = zeros(N_pointsL*N_pointsP*N_heights,1); %// viewing angle
THETA_INC = zeros(N_pointsL*N_pointsP*N_heights,1); %// inc. angle
AMATRIX   = zeros(N_pointsL*N_pointsP*N_heights,N_coeffs); %// design matrix

for k = 0:(N_heights-1)
    HEIGHT = p_baseline.H_min + k*deltaheight;
    ELLIPS.a = ell.a + HEIGHT;
    ELLIPS.b = ell.b + HEIGHT;
    if nargin < 3
        for i = 0:(N_pointsL-1) %azimuth lines
            line = current_window.top + i * deltalines;
            %______ Azimuth time for this line ______
            m_tazi = line2ta(line, m_info.sensor_para.t_Azi1, m_info.sensor_para.line_time_interval);
            
            %______ xyz for master satellite from time ______
            M = getxyz(m_orbit, m_tazi);
            
            %______ Loop over a pixels to compute baseline param. ______
            for j = 0:(N_pointsP-1) %range pixels
                
                %______ Range time for this pixel ______
                
                pixel = current_window.left + j * deltapixels;
                P = lph2xyz(m_orbit, m_info.sensor_para, line, pixel, 0);
                
                %______ Compute xyz for slave satellite from P ______
                slave= xyz2t(P, s_info,s_orbit);
                
                %_________________ Slave position ____________________
                S = getxyz(s_orbit, slave.timeAzimuth);
                
                %______ Compute angle between near parallel orbits ______
                %Mdot = getxyzdot(m_orbit, m_tazi);
                %Sdot = getxyzdot(s_orbit, s_tazi);
                %angleorbits = vec_angle(Mdot, Sdot);
                %p_baseline.orbit_convergence = angleorbits;
%               disp('Angle between orbits master-slave (at l,p) =');
%               disp([num2str(line) num2str(pixel) '=' num2str(rad2deg(angleorbits))]);
                
                [B, Bpar, Bperp, theta] = BBparBperpTheta(M, P, S);
                theta_inc = IncidenceAngle(M, P);
                
                BPERP(cnt, 1) = Bperp;
                BPAR(cnt, 1) = Bpar;
                THETA(cnt, 1) = theta;
                THETA_INC(cnt, 1) = theta_inc;
                
%               // --- B(l,p,h) = a000 + 
%               //                a100*l   + a010*p   + a001*h   +
%               //                a110*l*p + a101*l*h + a011*p*h +
%               //                a200*l^2 + a020*p^2 + a002*h^2
                
                AMATRIX(cnt,1)   =  1.0;
                AMATRIX(cnt,2)   = normalize(line,p_baseline.L_min,p_baseline.L_max);
                AMATRIX(cnt,3)   = normalize(pixel,p_baseline.P_min,p_baseline.P_max);
                AMATRIX(cnt,4)   = normalize(HEIGHT,p_baseline.H_min,p_baseline.H_max);
                AMATRIX(cnt,5)   = normalize(line,p_baseline.L_min,p_baseline.L_max) * ...
                    normalize(pixel,p_baseline.P_min,p_baseline.P_max);
                AMATRIX(cnt,6)   = normalize(line,p_baseline.L_min,p_baseline.L_max) * ...
                    normalize(HEIGHT,p_baseline.H_min,p_baseline.H_max);
                AMATRIX(cnt,7)   = normalize(pixel,p_baseline.P_min,p_baseline.P_max)* ...
                    normalize(HEIGHT,p_baseline.H_min,p_baseline.H_max);
                AMATRIX(cnt,8)   = (normalize(line,p_baseline.L_min,p_baseline.L_max))^2;
                AMATRIX(cnt,9)   = (normalize(pixel,p_baseline.P_min,p_baseline.P_max))^2;
                AMATRIX(cnt,10)  = (normalize(HEIGHT,p_baseline.H_min,p_baseline.H_max))^2;
                cnt = cnt + 1;
            end%j = 0:(N_pointsP-1) %range pixels
        end%for i = 0:(N_pointsL-1) %azimuth lines
    else%vargin = 3
        disp('depend on coregistration model, not finish yet!');
    end%vargin < 3
end%k = 0:(N_heights-1)
    
    N = AMATRIX' * AMATRIX;
    rhsBperp = AMATRIX' * BPERP;
    rhsBpar  = AMATRIX' * BPAR;
    rhsT     = AMATRIX' * THETA;
    rhsT_INC = AMATRIX' * THETA_INC;
%     Qx_hat = inv(N);
    
    rhsBperp = N \ rhsBperp;
    rhsBpar  = N \ rhsBpar;
    rhsT     = N \ rhsT;
    rhsT_INC = N \ rhsT_INC;
    
    maxdev = max(max(abs(N \ N - eye(length(N)))));
         
     if maxdev > 0.01
         error(['baseline: maximum deviation N*inv(N) from unity = ' num2str(maxdev)]);
     elseif maxdev > 0.001
         warning(['baseline: maximum deviation N*inv(N) from unity = ' num2str(maxdev)]);
     end
     % === Copy estimated coefficients to private members ===
     p_baseline.BPERP_cf     = rhsBperp;
     p_baseline.BPAR_cf      = rhsBpar;
     p_baseline.THETA_cf     = rhsT;
     p_baseline.THETA_INC_cf = rhsT_INC;
     %p_baseline = class(p_baseline,'baseline');
end   
    
    
%    /****************************************************************
%    * Return baselineparameters                                    *
%    ****************************************************************/

function [B, Bpar, Bperp, theta] = BBparBperpTheta(M, P, S)
Master=[M.X,M.Y,M.Z];
Point=[P.X,P.Y,P.Z];
Slave=[S.X,S.Y,S.Z];
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
% ****************************************************************
%  * returns incidence angle in radians based on coordinate of    *
%  * point P on ellips and point M in orbit                       *
%  #%// Bert Kampes, 06-Apr-2005
%  ****************************************************************/
        
    function inc_angle = IncidenceAngle(M, P)
    Master(1)=M.X;Master(2)=M.Y;Master(3)=M.Z;
Point(1)=P.X;Point(2)=P.Y;Point(3)=P.Z;
        r1 = Master - Point;
        inc_angle = vec_angle(Point, r1);
    end        
  
    %angle of two vectors
    function angle_of_vec = vec_angle(vec1, vec2)
    %vex1(1)=vec1.X;vex1(2)=vec1.Y;vex1(3)=vec1(3);
    %vex2(1)=vec2.X;vex2(2)=vec2.Y;vex2(3)=vec2(3);
        in_vec = vec1 * vec2';
        angle_of_vec = acos(in_vec / (norm(vec1)*norm(vec2)));
    end
function [POS]=getxyz(coeff,azTime)
azTimeNormal = (azTime - coeff.time(floor(length(coeff.time) / 2)+1))./10;
%azTimeNormal = azTime;
POS.X=polyVal1D(azTimeNormal, coeff.X);
POS.Y=polyVal1D(azTimeNormal, coeff.Y);
POS.Z=polyVal1D(azTimeNormal, coeff.Z);
end
function [sum]=polyVal1D( aziTime,coeffs)
sum = 0.0;
d=length(coeffs);
        while d > 0
            sum =sum.* aziTime;
            sum =sum+ coeffs(d);
            d=d-1;
        end
end
function [POS]=getxyzdot(coeff,azTime)
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
            
        
    
    
        
    
        
    
    








 


