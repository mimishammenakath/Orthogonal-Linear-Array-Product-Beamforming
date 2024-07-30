%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vertical and Horizontal Linear Array Beamforming using a cross array                      %                               
%  Parameters:                                                                                                                                %
%  1. Sensor_data_V : Sensor data of Vertical linear array
%  2. Sensor_data_H : Sensor data of Horizontal linear array
%  3. timeVector    : Time vector corresponding to the data
%  4. delay_V       : Delay for the Vertical linear array
%  5. delay_H       : Delay for the Horizontal linear array
%  6. range_index   : Time index corresponding to the range
%  7. beamformed_data : Beamformed output for a 2D slice at the
%  corresponding range
% Author: Mimisha M Menakath                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beamformed_data=VHLABF_Cross(Sensor_data_V,Sensor_data_H,timeVector,delay_V,delay_H,range_index)

InputParameters;

%% Proposed beamforming
beamformed_data_V=LABF_Cross(Sensor_data_V,timeVector,delay_V,range_index,N,Mb,Nb);
beamformed_data_H=LABF_Cross(Sensor_data_H,timeVector,delay_H,range_index,M,Mb,Nb);
beamformed_data=beamformed_data_V.*beamformed_data_H;
beamformed_data=sqrt(abs(beamformed_data));
end