%% Input Paramerers
M=24;                                          %number of sensor elements along x direction
N=24;                                           %number of sensor elements along y direction

Mb=120;                                         %number of beams in azimuth direction
Nb=Mb;                                         %number of beams in elevation direction
c=1500;                                        %speed of sound in the medium in m/s
                                        
alpha_i=-30*pi/180;                            %inital azimuth angle in radian
alpha_f=-alpha_i;                              %final azimuth angle in radian

beta_i=-30*pi/180;                             %inital elevation angle in radian
beta_f=-beta_i;                                %final elevation angle in radian

 
ds_alpha=(alpha_f-alpha_i)/(Mb);             %beam spacing along azimuth direction
ds_beta=(beta_f-beta_i)/(Nb);                %beam spacing along elevation direction
s_alpha=(sin(alpha_f)-sin(alpha_i))/(Mb);    %beam spacing along azimuth direction
s_beta=(sin(beta_f)-sin(beta_i))/(Nb);       %beam spacing along elevation direction