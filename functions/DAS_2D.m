%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conventional delay and sum beamforming                      %                               
%  Parameters:                                                                                                                                %
%  1. rawData       : Sensor data of the planar array
%  2. timeVector    : Time vector corresponding to the data
%  3. delay         : Delay for the planar array
%  4. range_index   : Time index corresponding to the range
%  5. beamformed_data : Beamformed output for a 2D slice at the
%  corresponding range
% Author: Mimisha M Menakath and Mahesh Raveendranatha Panicker                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beamformed_data=DAS_2D(rawData,timeVector,delay,range_index)

InputParameters;
r=range_index;
%% DAS beamforming
delay_compensation=zeros(M*N,Mb*Nb);
delay_compensated_signal=zeros(M*N,Mb*Nb);
delay=reshape(delay,M*N,Mb*Nb);
for Nc=1:M*N
    delay_compensation(Nc,:)=timeVector(r)-delay(Nc,:);
    delay_compensated_signal(Nc,:)=interp1(timeVector,(rawData(:,Nc))',delay_compensation(Nc,:),'spline',0);
end

beamformed_data=sum(delay_compensated_signal);
beamformed_data=reshape(beamformed_data,Nb,Mb);
beamformed_data=beamformed_data';
end
