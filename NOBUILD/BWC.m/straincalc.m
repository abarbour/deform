%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%straincalc.m
%This takes in the parameters to be used in strain.m and outputs the strain
%data
%Written by Brendan Crowell - bwcrowel@ucsd.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variables
dt = 0.1;%grid spacing
aa = 7;%grid distance weighting value, in km
minlat = 32.8;%minimum latitude
maxlat = 33.5;%maximum latitude
minlon = -116;%minimum longitude
maxlon = -115.2;%maximum latitude
lats = minlat:dt:maxlat;
lons = minlon:dt:maxlon;
[LAT,LON]=meshgrid(lats,lons);
[lon,lat,e,n,ee,ne,~,site]=textread('vel_diff.txt','%f %f %f %f %f %f %f %s');%gps data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Strain
[S1,S2,SS,DIL,ROT,ANG,ERR]=strain(lon,lat,e,n,LON,LAT,aa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write out the strain
fid = fopen('strain_diff.txt','wt');
fprintf(fid,'%s %s %s %s %s %s %s %s %s %s\n',...
    '#','lon','lat','S1','S2','MaxShearStrain','Dilatation','Rotation','AngleMaxShear','Errors');
for i=1:length(LON(:,1))
    for j = 1:length(LON(1,:))
        fprintf(fid,'%1.2f %1.2f %1.5e %1.5e %1.5e %1.5e %1.5e %1.2f %1.5e\n',...
            LON(i,j),LAT(i,j),S1(i,j),S2(i,j),(S1(i,j)-S2(i,j))/2,DIL(i,j),ROT(i,j),ANG(i,j)*180/pi,ERR(i,j));
    end
end
fclose(fid);




