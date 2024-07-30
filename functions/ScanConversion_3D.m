%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Scan Conversion                   %                               
%  Parameters:                                                                                                                                %
%  1. alpha_vec: Vector of Scaned azimuth angles
%  2. beta_vec : Vector of Scaned elevation angles
%  3. r_vec    : Vector of ranges 
%  4. image_3d : 3D image in polar cocordinate system
%  5. x_res    : X resolution
%  6. y_res    : Y resolution
%  7. z_res    : Z resolution
%  8. b_mode   : Scanconverted b_mode image in cartesian coordinate system
% Author: Mimisha M Menakath                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scanconversion
function b_mode= ScanConversion_3D(alpha_vec,beta_vec,r_vec,image_3d,x_res,y_res,z_res)

[alpha,beta,r]=ndgrid(alpha_vec,beta_vec,r_vec);

x_grid=r.*sin(alpha);
y_grid=r.*sin(beta);
z_grid=r.*sqrt((cos(alpha)).^2-(sin(beta)).^2);

x_vec=linspace(min(x_grid(:)),max(x_grid(:)),x_res);
y_vec=linspace(min(y_grid(:)),max(y_grid(:)),y_res);
z_vec=linspace(min(z_grid(:)),max(z_grid(:)),z_res);

[x_q,y_q,z_q]=ndgrid(x_vec,y_vec,z_vec);


%cartetian to polar
r_new=sqrt(x_q.^2+y_q.^2+z_q.^2);
alpha_new=asin(x_q./r_new);
beta_new=asin(y_q./r_new);
 
b_mode = interp3(alpha_vec ,beta_vec, r_vec ,image_3d,alpha_new,beta_new ,r_new, 'linear',0);
b_mode(isinf(b_mode))=0;
b_mode(isnan(b_mode))=0;
end


