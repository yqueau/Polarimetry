%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of Mueller imaging polarimetry after calibration based on 
%   [Quéau et al., 'Design and simplified calibration of a Mueller imaging
%  polarimeter for material classification', Opt. Lett. 43(20), 2018]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters which should be changed

% Path to the recorded images
images_path = 'Data/Smileys'; % Folder containing the images. Subfolder should be like this: $$lambda$$ / $$alpha_g$$_$$alpha_a$$.$$image_format$$
images_format = 'png'; % Format of the images
border_size = 0; % We will remove all intensity measurements on a border of this size

% Path to the calibration file
calibration_file = 'calibration_reflexion.mat';

% Controlled angles of the PSG and PSA retarders 
tab_alpha_g = 0:22.5:157.5; % Azimuts of the PSG retarder used for calibration
tab_alpha_a = 0:22.5:157.5; % Azimuts of the PSA retarder used for calibration

% Working wavelength
lambda = 630;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the images and store in a huge tensor I_all, 
path = [images_path '/' num2str(lambda)];
I_all = []; % Will be the stack of images: nalpha_g x nalpha_a x nrows x ncols 
I_all_vec = []; % Will be vectorized: nalpha_g.nalpha_a x nrows.ncols
ind = 0;
for i=1:length(tab_alpha_g)
	for j=1:length(tab_alpha_a)
		ind = ind+1;
		nom=[path '/' num2str(round(tab_alpha_g(i))) '_' num2str(round(tab_alpha_a(j))) '.' images_format];
		try
			mat = imread(nom);
			mat = mat(1+border_size:end-border_size,1+border_size:end-border_size); % Remove dark borders
		catch
			disp('Error: Some image non existing - replacing by NaN...');
			mat = NaN;
		end
		I_all(j,i,:,:) = mat;
		I_all_vec(ind,:) = mat(:);
	end
end
[~,~,nrows,ncols] = size(I_all);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate Mueller matrix field
load(calibration_file);

% Get delays of the PSG and PSA at the working wavelength
delta_g = k_1_g/lambda + k_2_g/(lambda^3); % Delay of the PSG
delta_a = k_1_a/lambda + k_2_a/(lambda^3); % Delay of the PSA

% Construct the PSA and PSG matrices such that I = A M G, in all pixels with I 8x8, A 8x4, M 4x4 and G 4x8
A = [ones(length(tab_alpha_a),1),...
            1-(1-cos(delta_a))*(sind(2*tab_alpha_a')).^2,...
            (1-cos(delta_a))*(cosd(2*tab_alpha_a')).*(sind(2*tab_alpha_a')),...
            sin(delta_a)*(sind(2*tab_alpha_a'))]; % Analyzer matrix
G = [ones(1,length(tab_alpha_g));...
	1-(1-cos(delta_g))*(sind(2*tab_alpha_g)).^2;...
	(1-cos(delta_g))*(cosd(2*tab_alpha_g)).*(sind(2*tab_alpha_g));...
	-sin(delta_g)*(sind(2*tab_alpha_g))]; % Generator matrix

% Construct the vectorized form C such that I = C M, in all pixels with I 64x1, C 64x16 and M 16x1        
C = kron(G',A);        

% Solve the overdetermined linear system
M_field = C\I_all_vec; 

% Turn into a matrix field
M_field = reshape(M_field',nrows,ncols,4,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the estimated field
figure('units','normalized','outerposition',[0 0 1 1])
cpt = 0;
for i = 1:4
	for j = 1:4
	
		cpt = cpt+1;
		subplot(4,4,cpt)
	
		imagesc(M_field(:,:,i,j)./M_field(:,:,1,1))
		axis image
		axis tight
		axis off
		colorbar
		name = sprintf('$$M_{%d%d}/M_{11}$$',i,j);
		title(name,'Interpreter','Latex')
		set(gca,'Fontsize',24)
		
	end
end
drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


