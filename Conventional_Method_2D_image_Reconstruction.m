%-- Script for reconstructing one slice of a 3D image at a given range
%-- using conventional delay and sum beamforming
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

%% Beamforming
beamformed_data=zeros(Mb,Nb);
tic
beamformed_data=DAS_2D(rawData,timeVector,delay,range_index);
toc
beamformed_data=beamformed_data./(M*N);
Normalized_beamformed_data=abs(beamformed_data)/max(abs(beamformed_data(:)));

alpha=(alpha_i:ds_alpha:alpha_f-ds_alpha)*180/pi;
beta=(beta_i:ds_beta:beta_f-ds_alpha)*180/pi;
figure,imagesc(alpha,beta,Normalized_beamformed_data); ylabel('Elevation Angle (Deg)','FontSize',14,'FontWeight','bold') ;
xlabel('Azimuth Angle (Deg)','FontSize',35,'FontWeight','bold') ;colorbar;%caxis ([-0.1 2]);
set(gca,'FontSize',14);axis equal; axis tight;
