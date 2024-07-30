%-- Script for reconstructing one slice of a 3D image at a given range
%-- using orthogonal linear array product beamforming
%-- Authors: Mimisha M Menakath and Mahesh Raveendranatha Panicker
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 29 - July - 2024
%-------------------------------------------------------------------------%

clearvars;
clc;
close all;

addpath('functions');
addpath('Data')
%% input parameters

load('TestData.mat');
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
delay_V11=squeeze( delay(1,:,1:ceil(Mb/2),1:ceil(Nb/2)));  % delay for V1 & the 1st quadrant
delay_H11=squeeze( delay(:,1,1:ceil(Mb/2),1:ceil(Nb/2)));  % delay for H1 & the 1st quadrant
delay_V12= squeeze(delay(1,:,1:ceil(Mb/2),ceil(Nb/2)+1:Nb)); % delay for V1 & the 2nd quadrant
delay_H22= squeeze(delay(:,end,1:ceil(Mb/2),ceil(Nb/2)+1:Nb));% delay for H2 & the 2nd quadrant
delay_V23= squeeze(delay(end,:,ceil(Mb/2)+1:Mb,1:ceil(Nb/2)));% delay for V2 & the 3rd quadrant
delay_H13= squeeze(delay(:,1,ceil(Mb/2)+1:Mb,1:ceil(Nb/2)));  % delay for H1 & the 3rd quadrant
delay_V24= squeeze(delay(end,:,ceil(Mb/2)+1:Mb,ceil(Nb/2)+1:Nb)); % delay for V2 & the 4th quadrant
delay_H24= squeeze(delay(:,end,ceil(Mb/2)+1:Mb,ceil(Nb/2)+1:Nb)); % delay for H2 & the 4th quadrant



Sensor_data_V1=squeeze(Sensor_data(:,1,:)); % Sensor data for V1
Sensor_data_H1=squeeze(Sensor_data(:,:,1)); % Sensor data for H1
Sensor_data_V2=squeeze(Sensor_data(:,end,:)); % Sensor data for V2
Sensor_data_H2=squeeze(Sensor_data(:,:,end)); % Sensor data for H2
%% Beamforming
beamformed_data=zeros(Mb,Nb);
tic
beamformed_data(1:ceil(Mb/2),1:ceil(Nb/2))=VHLABF(Sensor_data_V1,Sensor_data_H1,timeVector,delay_V11,delay_H11,range_index);
toc
beamformed_data(ceil(Mb/2)+1:Mb,1:ceil(Nb/2))=VHLABF(Sensor_data_V1,Sensor_data_H2,timeVector,delay_V12,delay_H22,range_index);
beamformed_data(1:ceil(Mb/2),ceil(Nb/2)+1:Nb)=VHLABF(Sensor_data_V2,Sensor_data_H1,timeVector,delay_V23,delay_H13,range_index);
beamformed_data(ceil(Mb/2)+1:Mb,ceil(Nb/2)+1:Nb)=VHLABF(Sensor_data_V2,Sensor_data_H2,timeVector,delay_V24,delay_H24,range_index);

beamformed_data=beamformed_data./(M*N);
Normalized_beamformed_data=abs(beamformed_data)/max(abs(beamformed_data(:)));

alpha=(alpha_i:ds_alpha:alpha_f-ds_alpha)*180/pi;
beta=(beta_i:ds_beta:beta_f-ds_alpha)*180/pi;
figure,imagesc(alpha,beta,Normalized_beamformed_data); ylabel('Elevation Angle (Deg)','FontSize',14,'FontWeight','bold') ;
xlabel('Azimuth Angle (Deg)','FontSize',35,'FontWeight','bold') ;colorbar;%caxis ([-0.1 2]);
set(gca,'FontSize',14);axis equal; axis tight;
