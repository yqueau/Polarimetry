function tab_I = load_images(images_path,images_format,border_size,tab_theta_a,tab_alpha_g,tab_alpha_a,tab_lambda); 

	nb_theta_a = length(tab_theta_a); % Number of measurements in PSA polarizer
	nb_alpha_g = length(tab_alpha_g); % Number of measurements in PSG retarder
	nb_alpha_a = length(tab_alpha_a); % Number of measurements in PSA retarder
	nb_wavelengths = length(tab_lambda); % Number of wavelengths

	nb_images = nb_theta_a*nb_alpha_g*nb_alpha_a*nb_wavelengths; % Total number of images
	
	% Allocate memory for the measurements
	tab_I = NaN*ones(1,nb_images);

	% Load data
	ind = 0;
	progressbar('Loading images')
	for i=1:nb_wavelengths
		for j=1:nb_theta_a
			for k=1:nb_alpha_g
				for l=1:nb_alpha_a
					ind = ind+1;
					progressbar(ind/nb_images)
					file = sprintf('%s/%d/%d_%d_%d.%s',images_path,tab_lambda(i),round(rad2deg(tab_theta_a(j))),round(rad2deg(tab_alpha_g(k))),round(rad2deg(tab_alpha_a(l))),images_format);
					try
						I_im = imread(file);
						I_im = double(I_im(border_size:end-border_size,border_size:end-border_size)); % Remove dark borders
						val_im = median(I_im(:)); % Average values inside
						if(val_im==255)
							val_im = NaN
						end
					catch
						val_im = NaN
					end
					tab_I(ind) = val_im;
				end
			end
		end
	end
	progressbar(1)

end
