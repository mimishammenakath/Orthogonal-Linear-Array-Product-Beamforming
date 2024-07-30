%-------------------------------------------------------------------------%
%-- Script for simulating planar data of a point source in farfield for
%-- underwater acoustical 3D imaging
%-- Simulation Method: Analytical modelling.
%-- Authors: Mimisha M Menakath and Hareesh G
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 29 - July - 2024
%-------------------------------------------------------------------------%

clear
clc
close all
addpath('functions');
%% Planar array point source simulation for an Active sonar

% Speed of sound in water
c=1500; 

% Design Frequency
fc=500e3;

% Wavelength
lambda=c/fc;

% Inter element spacing
dR=0.499*lambda;

% Number of horizontal elements
NH=24;

% Number of Vertical elements

NV=24;

% sampling Frequency
Fs = 4*fc;

% Signal to Noise Ratio
SNR = 100;


%  Array shading (Apodization)
win=rectwin(NH*NV)./(NH*NV);

%% Simulation parameters 
% Frequency of Operation

% Lower Frequency
Fl = 480e3;

% Upper Frequency
Fh = 520e3;

% Pulse Width
pw=10e-3;

% Pulse Repition PRI
pri=0.15; % sec

FreqInfo = [Fl Fh Fs c pw,pri];

% Total signal Length
SigLen=pri*Fs;

% Lower Frequency for processing
Bbf1 = Fl; 

% Upper Frequency for processing
Bbf2 = Fh; 
%%  Array and Target Parameter For simulating data

% Sensor position simulation
[SensorPos]=UPA(NH,NV,dR);

% Target (Source) parameters
targetAzimuth =0;  %  target azimuth angle in degree
targetElevation =0; % target elevation angle in degree
targetRange=30; % m

%% Array Data Simulation
[signal,arrayData,delay] = SimulateData(targetElevation,targetAzimuth,targetRange,SNR,NH,NV,FreqInfo,dR,SensorPos);
TimeVector=((1:size(arrayData,1))-1).*1/Fs;
save('TestData','arrayData','fc','Fs','targetRange','TimeVector','dR','targetRange','NH','NV','SensorPos','-v7.3');
