function [senPosR]=UPA(NH,NV,dR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to get sensor position of a uniform planar array              %
%
% Developed by:
%                Mimisha M Menakath and Hareesh G

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=NH;% number of elements in x direction
N=NV;% number of elements in y direction
d=dR;% element spacing
zPos1 = zeros(M,N);          % Z position
xPos1 =([0:M-1]-(M-1)/2 )*d; % x position
yPos1=([0:N-1]-(N-1)/2 )*d;  % y position
[xPos,yPos]=meshgrid(xPos1,yPos1);
senPosR(:,1) = xPos(:); % x indices of all the sensor elemnts
senPosR(:,2) = yPos(:); % y indices of all the sensor elemnts
senPosR(:,3) = zPos1(:);% z indices of all the sensor elemnts

