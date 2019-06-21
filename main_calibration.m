%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bundle calibration of a polarimeter as described in
%   [Quéau et al., 'Design and simplified calibration of a Mueller imaging
%  polarimeter for material classification', Opt. Lett. 43(20), 2018]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

% Get some useful functions
addpath('Toolbox/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters which should be changed

% Path to the recorded images
images_path = 'Data/Reflexion/'; % Folder containing the images. Subfolder should be like this: $$lambda$$ / $$theta_a$_$$alpha_g$$_$$alpha_a$$.$$image_format$$
images_format = 'png'; % Format of the images
border_size = 150; % We will remove all intensity measurements on a border of this size, then take the median of the remaining intensities

% Set mode (reflexion or transmission)
reflexion = 1; % Set to 0 for transmission, 1 for reflexion

% Save results or not
do_save = 1; 
output_file = 'calibration_reflexion.mat';

% Controlled angles of the PSG and PSA 
tab_theta_a = deg2rad(0:45:315); % Angles of the PSA polarizer used for calibration
tab_alpha_g = deg2rad(0:22.5:157.5); % Azimuts of the PSG retarder used for calibration
tab_alpha_a = deg2rad(0:22.5:157.5); % Azimuts of the PSA retarder used for calibration
tab_lambda = [460 510 540 580 630 650]; % Wavelengths used for calibration

% Rough initial guess of the parameters to estimate
theta_0_a = deg2rad(45); % "Zero" of the PSA polarizer
alpha_0_g = deg2rad(45); % "Zero" of the PSG retarder
alpha_0_a = deg2rad(45); % "Zero" of the PSA retarder
k_1_g = 1000; % Cauchy's first parameter of the PSG delay
k_2_g = 1e7; % Cauchy's second parameter of the PSG delay
k_1_a = 1000; % Cauchy's first parameter of the PSG delay
k_2_a = 1e7; % Cauchy's second parameter of the PSG delay
b_i = 150*ones(1, length(tab_lambda)); % Scale parameters (sensor response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters which should not need to be changed

% Set optimization parameters;
options.lsqnonlin = optimoptions('lsqnonlin');
options.lsqnonlin.Algorithm = 'levenberg-marquardt';
options.lsqnonlin.Display = 'iter-detailed' ;
options.lsqnonlin.FunctionTolerance = 1e-9 ;
options.lsqnonlin.MaxIterations = 200 ;
options.lsqnonlin.OptimalityTolerance = 1e-9 ;
options.lsqnonlin.StepTolerance = 1e-9 ;
options.lsqnonlin.SpecifyObjectiveGradient = true ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration of the polarimeter

% Load a vector containing all intensity measurements
tab_I = load_images(images_path,images_format,border_size,tab_theta_a,tab_alpha_g,tab_alpha_a,tab_lambda); 

% Do bundle adjustment
x_in = [alpha_0_g k_1_g k_2_g alpha_0_a k_1_a k_2_a theta_0_a b_i]; % Initial guess
if(reflexion)
	[x_out] = lsqnonlin(@(x)calibration_bundle_mirror(x,tab_I,tab_lambda,tab_theta_a,tab_alpha_g,tab_alpha_a),x_in,[],[],options.lsqnonlin);		
else
	[x_out] = lsqnonlin(@(x)calibration_bundle(x,tab_I,tab_lambda,tab_theta_a,tab_alpha_g,tab_alpha_a),x_in,[],[],options.lsqnonlin);		
end

% Retrieve estimated zeros
alpha_0_g = rad2deg(mod(x_out(1),pi));
alpha_0_a = rad2deg(mod(x_out(4),pi));
theta_0_a = rad2deg(mod(x_out(7),pi));
% Retrieve estimated delay functions
k_1_g = x_out(2);
k_2_g = x_out(3);
k_1_a = x_out(5);
k_2_a = x_out(6);
% Retrieve estimated proportionality constants (probalby not interesting)
b_i = x_out(8:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display estimation results

disp(sprintf('Estimated theta_0_a: %.3f deg',theta_0_a));
disp(sprintf('Estimated alpha_0_g: %.3f deg',alpha_0_g));
disp(sprintf('Estimated alpha_0_a: %.3f deg',alpha_0_a));

% Draw the delay curves for verification
tab_lambda_curve = (min(tab_lambda)-10):(max(tab_lambda)+10);
delta_g = [1./tab_lambda' 1./(tab_lambda').^3]*[k_1_g; k_2_g]; % Measurements
tab_delta_g_curve = k_1_g./tab_lambda_curve+k_2_g./tab_lambda_curve.^3; % Curve fit
figure(1)
plot(tab_lambda,rad2deg(delta_g),'+b','Markersize',10)
hold on;
plot(tab_lambda_curve,rad2deg(tab_delta_g_curve),'-b')
xlabel('$$\lambda$$ (nm)','Interpreter','Latex','FontSize',20)
ylabel('$$\delta^{G}$$ (deg)','Interpreter','Latex','FontSize',20)
title('Estimated PSG delay function')
axis tight
hold off;

delta_a = [1./tab_lambda' 1./(tab_lambda').^3]*[k_1_a; k_2_a];
tab_delta_a_curve = k_1_a./tab_lambda_curve+k_2_a./tab_lambda_curve.^3;
figure(2)
plot(tab_lambda,rad2deg(delta_a),'+b','Markersize',10)
hold on;
plot(tab_lambda_curve,rad2deg(tab_delta_a_curve),'-b')
xlabel('$$\lambda$$ (nm)','Interpreter','Latex','FontSize',20)
ylabel('$$\delta^{A}$$ (deg)','Interpreter','Latex','FontSize',20)
title('Estimated PSA delay function')
axis tight
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
if(do_save)
	save(output_file,'theta_0_a','alpha_0_g','alpha_0_a','k_1_g','k_2_g','k_1_a','k_2_a','b_i');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



