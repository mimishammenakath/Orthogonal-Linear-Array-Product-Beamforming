%-- Script for modelling an Square shape underwater target in k-Wave 3D grid.
%-- Authors:Mimisha M Menakath and Mahesh Raveendranatha Panicker
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 29 - July - 2024
%-------------------------------------------------------------------------%
%-- Dependencies:
%-- k-Wave Toolbox(http://www.k-wave.org/)
%-------------------------------------------------------------------------%
%-- Acknowledgements for smoothn function :
%-- 1. Garcia D. Robust smoothing of gridded data in one and higher dimensions
%-- with missing values. Computational Statistics & Data Analysis 2010; 54:1167-1178.
%-- 2.Garcia D. A fast all-in-one method for automated post-processing of PIV data.
%-- Exp Fluids, 2011;50:1247-1259. http://www.biomecardio.com/publis/expfluids11.pdf
%-------------------------------------------------------------------------%

clearvars;
clear all;
close all;
% =========================================================================
% DEFINE LITERALS
% =========================================================================

% select which k-Wave code to run
%   1: MATLAB CPU code
%   2: MATLAB GPU code
%   3: C++ code
%   4: CUDA code
model           = 4;

% medium parameters
c0              = 1500;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]

% source parameters
source_f0       = 500e3;      % source frequency [Hz]
source_roc      = 20e-3;    % bowl radius of curvature [m]
source_diameter =20e-3;    % bowl aperture diameter [m]
source_amp      = 1e6;      % source pressure [Pa]

% grid parameters
axial_size      = 600e-3;    % total grid size in the axial dimension [m]
lateral_size    =250e-3;    % total grid size in the lateral dimension [m]
% calculate the grid spacing based on the PPW and F0
ppw             = 3;        % number of points per wavelength
dx = c0 / (ppw * source_f0);
dy=dx;
dz=dx;% [m]
% compute the size of the grid
source_x_offset = 20;
Nx = roundEven(axial_size / dx) + source_x_offset;
Ny = roundEven(lateral_size / dx);
Nz = Ny;

bli_tolerance   = 0.01;     % tolerance for truncation of the off-grid source points
upsampling_rate = 3;       % density of integration points relative to grid

% =========================================================================
% RUN SIMULATION
% =========================================================================

% --------------------
% GRID
% --------------------


% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);
c0=1500;                                                % Assumed average speed of sound
cfl                 = [];
t_end = (Nx * dx) * 2.2/ c0;

kgrid.makeTime(c0, [], t_end);

% --------------------
% SOURCE
% --------------------
source_cycles       = 3;        % number of tone burst cycles
source_sig  = source_amp*toneBurst(1/kgrid.dt, source_f0 , source_cycles);

% set bowl position and orientation
bowl_pos = [kgrid.x_vec(1) + source_x_offset * kgrid.dx, 0, 0];
% focus_pos = [kgrid.x_vec(end), 0, 0];
focus_pos = [kgrid.x_vec(1), 0,0];
% create empty kWaveArray
karray = kWaveArray('BLITolerance', bli_tolerance, 'UpsamplingRate', upsampling_rate, 'SinglePrecision', true);

% add bowl shaped element
karray.addBowlElement(bowl_pos, source_roc, source_diameter, focus_pos);

% assign binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% assign source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

% --------------------
% MEDIUM
% --------------------
rho0 = 1000;                    % [kg/m^3]
alpha_coeff = 0.0000;      % [dB/(MHz^y cm)]
alpha_power = 1.05;
%medium.BonA = 6;
background_map_mean = 1;
background_map_std = 0.001;
background_map = background_map_mean + background_map_std * randn([Nx, Ny, Nz]);

% assign medium properties
medium.sound_speed = c0;
medium.density = rho0;
medium.sound_speed            = c0*ones(Nx, Ny, Nz).*background_map;
medium.density                = rho0*ones(Nx, Ny, Nz).*background_map;

scattering_region4= zeros(Nx,Ny,Nz);
scattering_region4(500:550,20:40,20:120)=1;
scattering_region4(500:550,25:35,55:95)=0;

scattering_region3= zeros(Nx,Ny,Nz);
scattering_region3(500:550,20:120,20:40)=1;
scattering_region3(500:550,55:195,25:35)=0;

scattering_region2= zeros(Nx,Ny,Nz);
scattering_region2(500:550,20:120, 100:120)=1;
scattering_region2(500:550,55:95,105:115)=0;

scattering_region1= zeros(Nx,Ny,Nz);
scattering_region1(500:550,100:120, 20:120)=1;
scattering_region1(500:550,105:115,55:95)=0;



size_1=size(medium.density(scattering_region1 == 1) );
size_2=size(medium.density(scattering_region2 == 1) );
size_3=size(medium.density(scattering_region3 == 1) );
size_4=size(medium.density(scattering_region4 == 1) );

medium.sound_speed(scattering_region1 == 1) =2000*ones(size_1).*(1+0.05 *randn(size_1));%scattering_c0(scattering_region1 == 1);
medium.density(scattering_region1 == 1) = 3000*ones(size_1).*(1+0.05 *randn(size_1));
medium.sound_speed(scattering_region2 == 1) =2000*ones(size_2).*(1+0.05 *randn(size_2));%scattering_c0(scattering_region1 == 1);
medium.density(scattering_region2 == 1) = 3000*ones(size_2).*(1+0.05 *randn(size_2));
medium.sound_speed(scattering_region3 == 1) =2000*ones(size_3).*(1+0.05 *randn(size_3));%scattering_c0(scattering_region1 == 1);
medium.density(scattering_region3 == 1) = 3000*ones(size_3).*(1+0.05 *randn(size_3));
medium.sound_speed(scattering_region4 == 1) =2000*ones(size_4).*(1+0.05 *randn(size_4));%scattering_c0(scattering_region1 == 1);
medium.density(scattering_region4 == 1) = 3000*ones(size_4).*(1+0.05 *randn(size_4));


medium.sound_speed=smoothn(medium.sound_speed,0.1);
medium.density=smoothn(medium.density,0.1);


medium.alpha_coeff            = alpha_coeff*ones(Nx, Ny, Nz);
medium.alpha_power            = alpha_power;


% --------------------
% SENSOR
% --------------------
array_width=48; %transducer array width
array_length=48; %transducer array length
% set sensor mask to record central plane, not including the source point
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1,(Ny/2 - array_width/2):1:(Ny/2 - array_width/2)+47,(Nz/2)-array_length/2:1:(Nz/2)-array_length/2+47) = 1;

% record the pressure
sensor.record = {'p'};

% record only the final few periods when the field is in steady state
sensor.record_start_index=1;
% --------------------
% SIMULATION
% --------------------

% set input options
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};

% run code
switch model
    case 1

        % MATLAB CPU
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
            input_args{:}, ...
            'DataCast', 'single', ...
            'PlotScale', [-1, 1] * source_amp);

    case 2

        % MATLAB GPU
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
            input_args{:}, ...
            'DataCast', 'gpuArray-single', ...
            'PlotScale', [-1, 1] * source_amp);

    case 3

        % C++
        sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});

    case 4

        % C++/CUDA GPU
        sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});

end
sensor_data_filtered=filtfilt(fliplr(source_sig),1,double(sensor_data.p'));

x_grid=kgrid.y(:,:,:);
z_grid=kgrid.x(:,:,:);
y_grid=kgrid.z(:,:,:);
z_grid=z_grid-min(z_grid(:));
y_grid=y_grid-min(y_grid(:));
x_axis=linspace(min(kgrid.y(:)),max(kgrid.y(:)),size(x_grid,2));
y_axis=linspace(min(kgrid.z(:)),max(kgrid.z(:)),size(y_grid,3));
temp=x_axis(Ny/2-array_width/2:1:Ny/2-array_width/2+47);
width_probe_geometry=linspace(min(temp),max(temp),48);
temp=y_axis(Nz/2-array_length/2:1:Nz/2-array_length/2+47);
length_probe_geometry=linspace(min(temp),max(temp),48)';

timeVector=kgrid.t_array;
rawData=sensor_data_filtered';
fs=1/kgrid.dt;
d=dx;
fc=source_f0;
scattering_region=permute(single(scattering_region1|scattering_region2|scattering_region3|scattering_region4),[2 3 1]);
save('square','rawData','timeVector','fs','d','dx','fc','Nx','Ny','Nz','width_probe_geometry','length_probe_geometry','scattering_region');