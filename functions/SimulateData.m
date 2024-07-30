function [signal1, SensorData, Delay] = SimulateData(Theta, phi,Range,SNR,NH,NV,FreqInfo,dR,SensorPos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to simulate planar data for a point source in farfield               %
%
% Developed by:
%                Mimisha M Menakath and Hareesh G
% Parameters
% Theta - Azimuth Angle
% phi -  Elevation Angle
% SensorPos - Sensor Position
% SensorAngs - Angular positions of sensors
% NH - Number of horizontal elements
% NV - Number of vertical element
% NoOfElementsforAzimuthBeam - Number of elements for forming beams

% Uncomment the commented lines for adding ambient noise to the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fl = FreqInfo(1);
Fh = FreqInfo(2);
Fs = FreqInfo(3);
c  = FreqInfo(4);
pw = FreqInfo(5);
pri= FreqInfo(6);
P= SensorPos;

SigLen=pri*Fs;
SensorData = zeros(SigLen,NH*NV);


N=12;
B = fir1(N,4.0*[Fl Fh]/(Fs*2))'; % Filter Coefficients.

A=1;
U = [sind(phi);sind(Theta); sqrt((cosd(phi)).^2-(sind(Theta)).^2)]; % direction cosines
%% delay calculation

delay=zeros(NH,NV);

for m=1:NH                             % sensor position in x direction
    for n=1:NV                         % sensor poistion in y direction
        delay(m,n)=(((m-1)-(NH-1)/2 )*dR*sind(phi)+((n-1)-(NV-1)/2) *dR*sind(Theta))/c;
    end
end
delay=reshape(delay,1,NH*NV);
Delay=delay;
Len = pw*Fs;
signal1 = A*cos(2*pi*(Fl+Fh)/2*(0:Len-1)./Fs);

signal=[zeros(1,100) signal1 zeros(1,100)];

L=length(signal);

%% shifting in frequency domain
%
% signal=repmat(signal',1,NH*NV);
%  signal_f = fft(signal);
%  signal_f1 = signal_f(1:L/2+1,:);
%  f = 0:Fs/L:Fs/2;
%
% %  for l=1:length(f)
%      for j=1:NH*NV
%        signal_f2=signal_f1';
% %  signal_f3(j,l) = signal_f2(j,l).*exp(-1i*2*pi*f(l).*Delay(j));
%  signal_f3(j,:) = signal_f2(j,:).*exp(-1i*2*pi*f.*Delay(j));
%      end
% %  end
% signal_f4 = [signal_f3 conj(fliplr(signal_f3(:,2:end-1)))];
%
%  signal_shifted = ifft(signal_f4','symmetric');
%
%  RangeBin=floor(2*(Range)/c*Fs);
%
%  SensorData(RangeBin-100:RangeBin-100+L-1,:)=signal_shifted;
%% shifting in time domain

t=(0:L-1)./Fs;
signal=repmat(signal,NH*NV,1);
for j=1:NH*NV
    delayedTime=t+Delay(j);
    signal_shifted(j,:)=interp1(t,signal(j,:),delayedTime,'spline',0);
end
RangeBin=floor(2*(Range)/c*Fs);
SensorData(RangeBin-100:RangeBin-100+L-1,:)=signal_shifted';

