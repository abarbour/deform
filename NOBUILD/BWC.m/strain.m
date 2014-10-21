%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%strain.m
%This program takes in a grid of points and velocities/displacements and
%outputs the strain components at each grid point.  It uses methods from
%Shen et al. [1996] and Allmendinger et al. [2007].  The variables are as 
%follows:
%lon - longitudes of measurements
%lat - latitudes of measurements
%e - east displacement/velocities of measurements in mm/yr
%n - north displacement/velocities of measurements in mm/yr
%LON - longitude grid points.  Use meshgrid to form this with LAT
%LAT - latitude grid points
%aa - alpha value of grid distance weight matrix laid out by Allmendinger
%et al. [2007]
%Written by Brendan Crowell - bwcrowel@ucsd.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S1,S2,SS,DIL,ROT,ANG,ERR]=strain(lon,lat,e,n,LON,LAT,aa)

[llat,llon]=degreelen(LAT);%lengths of degree of lat and lon

for i = 1:length(LON(:,1))
    for j = 1:length(LON(1,:))
        u=[];G=[];
        W=zeros(2*length(lon),2*length(lon));
        for k=1:length(lon)   
            dx1 = (lon(k)-LON(i,j))*llon(i,j);
            dx2 = (lat(k)-LAT(i,j))*llat(i,j);
            d = (dx1^2+dx2^2)^0.5/1000;%distance to each grid point
            u=[u;e(k)/1000;n(k)/1000];%displacements or velocities
            G=[G;1 0 dx1 dx2 0 0;0 1 0 0 dx1 dx2];%green's function representation of translation plus velocity gradient tensor components
            W(2*(k-1)+1,2*(k-1)+1)=exp(-d^2/2/aa^2);%Weighting matrix, for east
            W(2*(k-1)+2,2*(k-1)+2)=exp(-d^2/2/aa^2);%Weighting matrix, for north
        end

        [S] = lsqlin(G'*W*G,G'*W*u,[],[]);%solve for strain

        ERR(i,j) = std((G*S-u)./length(u)) /0.1 /(llat(i,j)^2+llon(i,j)^2)^0.5; %error estimate in each grid point

        SXX(i,j) = S(3);%terms of the velocity gradient tensor
        SYY(i,j) = S(5);
        SXY(i,j) = S(4);
        SYX(i,j) = S(6);
        
        SS(i,j) = 0.5.*(SXY(i,j)+SYX(i,j));%shear strain
        DIL(i,j) = SXX(i,j)+SYY(i,j);%dilatation
        S1(i,j) = (SXX(i,j)+SYY(i,j))./2+(((SXX(i,j)-SYY(i,j))./2).^2+SS(i,j).^2).^0.5;%First Principal Component - Mohr's Circle definition
        S2(i,j) = (SXX(i,j)+SYY(i,j))./2-(((SXX(i,j)-SYY(i,j))./2).^2+SS(i,j).^2).^0.5;%Second Principal Component
        
        if (S2(i,j) > S1(i,j))%Change the order if S2 is larger
            S2n = S1(i,j);
            S1(i,j) = S2(i,j);
            S2(i,j) = S2n;
        end
        ROT(i,j) = -0.5.*(SXY(i,j)-SYX(i,j)).*180./pi; %rotation rate, degrees per year
        ANG(i,j) = 0.5.*atan(-(SXX(i,j)-SYY(i,j))/2/SS(i,j)); %angle of plane of maximum shear strain - Mohr's Circle   
    end
end

return
