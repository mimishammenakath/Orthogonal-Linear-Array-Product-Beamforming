%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Linear Array Delay and Sum Beamforming for 3D                     %                               
%  Parameters:                                                                                                                                %
%  1. Sensor_data   : Sensor data of thelinear array
%  2. timeVector    : Time vector corresponding to the data
%  3. delay         : Delay for the linear array
%  4. range_index   : Time index corresponding to the range
%  5. M             : Number of elements in the linear array
%  6. Mb            : Number of beams in azimuth direction
%  7. Nb            : Number of beams in elevation direction
%  8. beamformed_data : Beamformed output for a 2D slice at the corresponding range
% Author: Mimisha M Menakath                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beamformed_data=LABF_Cross(Sensor_data,timeVector,delay,range_index,M,Mb,Nb)

r=range_index;
%%  DAS beamforming

delay_compensation=zeros(M,Mb*Nb);
delay_compensated_signal=zeros(M,Mb*Nb);
delay=reshape(delay,M,Mb*Nb);

for Nc=1:M
    delay_compensation(Nc,:)=timeVector(r)-delay(Nc,:);
    delay_compensated_signal(Nc,:)=interp1(timeVector,(Sensor_data(:,Nc))',delay_compensation(Nc,:),'spline',0);
end

beamformed_data=sum(delay_compensated_signal);
beamformed_data=reshape(beamformed_data,Mb,Nb);
beamformed_data=beamformed_data';
end

