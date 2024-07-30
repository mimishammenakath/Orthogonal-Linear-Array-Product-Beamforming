%-- Script for reconstructing 3D image of underwater targets using
%-- Conventional delay and sum beamforming
%-- Authors: Mimisha M Menakath and Mahesh Raveendranatha Panicker
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 29 - July - 2024
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-- Disclaimer: Please note that the 'InputParameters' need to be changed
%-- according to the data. For the given data, please change the parameters M
%-- and N as 48.
%-------------------------------------------------------------------------%
clearvars;
clc;
close all;
addpath('functions');
addpath('Data')
%% input parameters
load('L_Shape');
rawData1=zeros(size(rawData));
rawData1(:,500:end)=rawData(:,500:end);
rawData=rawData1';
InputParameters;
timeVector=timeVector;
Sensor_data=reshape(rawData,[],M,N);

%% delay calculation
delay=zeros(M,N,Mb,Nb);
% r=Range;
for p=1:Mb                                     % azimuth angle selection
    for q=1:Nb                                 % elevation angle selection
        for m=1:M                              % sensor position in x direction
            for n=1:N                          % sensor poistion in y direction
                sin_alpha_p=sin(alpha_i+(p-1)*ds_alpha);
                sin_beta_q=sin(beta_i+(q-1)*ds_beta);
x_m=((m-1)-(M-1)/2)*d;
y_n=((n-1)-(N-1)/2)*d;
                delay(m,n,p,q)=((x_m*sin_alpha_p+y_n*sin_beta_q)/c);
                
            end
        end
    end
end

%% 
r_res=0.0005 ;                                                       % range resolution in meter                         
range_resolution=ceil(r_res*fs/(c*0.5));                               % range resolution in number of samples
r_d = 0.5*c.*timeVector(1:range_resolution:end);   % range vector
r_vec=1:range_resolution:length(timeVector);

%% Beamforming
beamformed_data=zeros(Mb,Nb,length(r_vec));
for r=1:length(r_vec)
    range_index=r_vec(r);
beamformed_data(:,:,r)=DAS_2D(rawData,timeVector,delay,range_index);
end
beamformed_data=beamformed_data./(M*N);
envelope=abs(beamformed_data);

%% Scanconversion
alpha_vec=alpha_i:ds_alpha:alpha_f-ds_alpha;
beta_vec=beta_i:ds_beta:beta_f-ds_beta;
x_res=250;                 % Resolution in x direction
y_res=250;                 % Resolution in y direction
z_res=250;                 % Resolution in z direction
b_mode= ScanConversion_3D(alpha_vec,beta_vec,r_vec,envelope,x_res,y_res,z_res);
b_mode_sc =(b_mode./mean(mean((b_mode))))*8;


 %% different slices of the 3D image 
dz=(r_d(end)-r_d(1))/z_res;
dx=(r_d(end)*sin(alpha_f)-r_d(end)*sin(alpha_i))/x_res;
dy=(r_d(end)*sin(beta_f)-r_d(end)*sin(beta_i))/y_res;
x=r_d(end)*sin(alpha_i):dx:r_d(end)*sin(alpha_f);
y=r_d(end)*sin(beta_i):dy:r_d(end)*sin(beta_f);
z=r_d(1):dz:r_d(end);

% yz projection
z_max=max(b_mode,[],1);
b_mode_sc_1 =(z_max./mean(mean((z_max))))*30;
figure,image(y,z,squeeze(b_mode_sc_1(1,:,:))'); % Plotting the speed of sound map
title('YZ Projection');xlabel('Y distance (m)','FontSize',12,'FontWeight','bold') ;
ylabel('Range (m)','FontSize',12,'FontWeight','bold') ; axis equal; axis tight;set(gca,'FontSize',14);axis equal; axis tight;
colormap(hot);colorbar;

% xz projection
z_max=max(b_mode,[],2);
b_mode_sc_1 =(z_max./mean(mean((z_max))))*30;
figure,image(x,z,squeeze(b_mode_sc_1(:,1,:))'); % Plotting the speed of sound map
title('XZ Projection');xlabel('X distance (m)','FontSize',12,'FontWeight','bold') ;
ylabel('Range (m)','FontSize',12,'FontWeight','bold') ; axis equal; axis tight;set(gca,'FontSize',14);axis equal; axis tight;
colormap(hot);colorbar;

% xy projection
z_max=max(b_mode,[],3);
b_mode_sc_1 =(z_max./mean(mean((z_max))))*60;
figure,image(x,y,b_mode_sc_1); % Plotting the speed of sound map
title('XY Projection');xlabel('X distance (m)','FontSize',12,'FontWeight','bold') ;
ylabel('Y distance (m)','FontSize',12,'FontWeight','bold') ;set(gca,'FontSize',14);
colormap(hot);colorbar;axis equal; axis tight;

