%-------------------------------------------------------------------------%
%-- Script for plotting point spread function of a double length cross
%-- array (DLCA) and an L array placed at the center (CLA) of the planar array for the
%-- product beamforming
%-- Authors: Mimisha M Menakath
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 29 - July - 2024
%-------------------------------------------------------------------------%
%-- Please change M and N in the InputParameters as 48.
%-------------------------------------------------------------------------%
clearvars;
clc;
close all;

%% input parameters
load('TestData_48_Elemnts.mat');
rawData=arrayData;
InputParameters;
timeVector=TimeVector;
d=dR;
fs=Fs;
Range=targetRange;
M=NH;
N=NV;
range_index=floor(2*(Range)/c*Fs);
Sensor_data=reshape(rawData,[],M,N);

%% delay calculation
delay=zeros(M,N,Mb,Nb);
for p=1:Mb                                     % azimuth angle selection
    for q=1:Nb                                 % elevation angle selection
        for m=1:M                              % sensor position in x direction
            for n=1:N                          % sensor poistion in y direction

                sin_alpha_p=sin(alpha_i+(p-1)*ds_alpha);
                sin_beta_q=sin(beta_i+(q-1)*ds_beta);
                delay(m,n,p,q)=(((m-1)-(M-1)/2)*d*sin_alpha_p+ ((n-1)-(N-1)/2)*d*sin_beta_q)/c;
                
            end
       
        end
            
    end
end
%% Uncomment the  following for the PSF of DLCA 
% delay_V=squeeze(delay(M/2,:,:,:));  % delay for Cross array
% delay_H=squeeze(delay(:,N/2,:,:)); 
% Sensor_data_V=squeeze(Sensor_data(:,M/2,:));   % data for Cross array
% Sensor_data_H=squeeze(Sensor_data(:,:,N/2));

%% Uncomment the  following for the PSF of CLA 
delay_V=squeeze(delay(M/2,1:N/2,:,:));  % delay for L array
delay_H=squeeze(delay(1:M/2,N/2,:,:)); 
Sensor_data_V=squeeze(Sensor_data(:,M/2,1:N/2));   % data for L array
Sensor_data_H=squeeze(Sensor_data(:,1:M/2,N/2));


%% Beamforming

beamformed_data=zeros(Mb,Nb);
tic
beamformed_data=VHLABF_Cross_Farfield(Sensor_data_V,Sensor_data_H,timeVector,delay_V,delay_H,range_index);
toc

beamformed_data=(beamformed_data)./(M*N);

Normalized_beamformed_data=abs(beamformed_data)/max(abs(beamformed_data(:)));
alpha=(alpha_i:ds_alpha:alpha_f-ds_alpha)*180/pi;
beta=(beta_i:ds_beta:beta_f-ds_alpha)*180/pi;

figure,imagesc(alpha,beta,Normalized_beamformed_data); ylabel('Elevation Angle (Deg)','FontSize',35,'FontWeight','bold') ;
xlabel('Azimuth Angle (Deg)','FontSize',14,'FontWeight','bold') ;colorbar;%caxis ([-0.1 2]);
set(gca,'FontSize',14);axis equal; axis tight;

indxAz=61;
indxEl=61;
%% PSF: Power Pattern
figure
subplot(2,1,1)
plot(alpha,20*log10(abs(beamformed_data(:,indxEl)./max(beamformed_data(:,indxEl)))));
xlabel('Azimuth Angle (Deg)','FontSize',35) ;
ylabel('Normalized power (dB)','FontSize',35) ;set(gca,'FontSize',30)
title("PSF: Azimuth Pattern")

%azimuth Beam width
[indx]=find(20*log10(abs(beamformed_data(:,indxEl)./max(beamformed_data(:,indxEl))))>-3);
Azimuth_BW=-2*(alpha(indx(1)))

subplot(2,1,2)
plot(beta,20*log10(abs(beamformed_data(indxAz,:)./max(beamformed_data(indxAz,:)))));
xlabel('Elevation Angle (Deg)','FontSize',35) ;
ylabel('Normalized power (dB)','FontSize', 35) ;set(gca,'FontSize',30)
title("PSF: Elevation Pattern")

%elevation Beam width
[indx]=find(20*log10(abs(beamformed_data(indxAz,:)./max(beamformed_data(indxAz,:))))>-3);
Elevation_BW=-2*(beta(indx(1)))