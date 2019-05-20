function Result=BLH2XYZ(B,L,H)  
H=double(H);
AE_WGS84=6378137.0;        
e2=0.00669438003551279091;   
N=AE_WGS84./sqrt(1.0-e2*sin(B/180*pi).*sin(B/180*pi)); 
TmpV=N+H; 
Result.X=TmpV.*cos(B/180*pi).*cos(L/180*pi); 
Result.Y=TmpV.*cos(B/180*pi).*sin(L/180*pi); 
Result.Z=(N*(1-e2)+H).*sin(B/180*pi); 
return