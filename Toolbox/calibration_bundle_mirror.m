function [F,J] = calibration_bundle_mirror(x, I, lambda, theta_a, alpha_g, alpha_a)
    
    alpha_0_g = x(1);
    k_1_g = x(2);
    k_2_g = x(3);
    alpha_0_a = x(4);
    k_1_a = x(5);
    k_2_a = x(6);
    theta_0_a = x(7);
    b = x(8:end);

    nb_wavelengths = length(lambda);
    nb_theta = length(theta_a);
    nb_alpha_g = length(alpha_g);
    nb_alpha_a = length(alpha_a);

    A = zeros(1,4);
    G = zeros(4,1);
    F = zeros(size(I));
    
    if (nargout > 1)
        J = zeros([length(I(:)) 1]);
        dG_dalpha0g = zeros(4,1);
        dG_dk1g = zeros(4,1);
        dG_dk2g = zeros(4,1);
        dA_dalpha0a = zeros(1,4);
        dA_dk1a = zeros(1,4);
        dA_dk2a = zeros(1,4);
        dA_dtheta0a = zeros(1,4);
    end
    
    ind = 0;
    for i=1:nb_wavelengths
        for j=1:nb_theta
            for k=1:nb_alpha_g
                for l=1:nb_alpha_a
                    ind = ind+1;
                    
                    if ~isnan(I(ind))
                        delta_g = k_1_g/lambda(i)+k_2_g/(lambda(i)^3);
                        delta_a = k_1_a/lambda(i)+k_2_a/(lambda(i)^3);
                        G(1) = 1;
                        G(2) = cos(delta_g)*sin(2*(alpha_g(k)-alpha_0_g))^2+cos(2*(alpha_g(k)-alpha_0_g))^2;
                        G(3) = (1-cos(delta_g))*cos(2*(alpha_g(k)-alpha_0_g))*sin(2*(alpha_g(k)-alpha_0_g));
                        G(4) = -sin(delta_g)*sin(2*(alpha_g(k)-alpha_0_g));
                        
                        A(1) = 1;
                        A(2) = cos(2*(theta_a(j)-theta_0_a))*(cos(2*(alpha_a(l)-alpha_0_a))^2+cos(delta_a)*sin(2*(alpha_a(l)-alpha_0_a))^2)+sin(2*(theta_a(j)-theta_0_a))*(1-cos(delta_a))*cos(2*(alpha_a(l)-alpha_0_a))*sin(2*(alpha_a(l)-alpha_0_a));
                        A(3) = cos(2*(theta_a(j)-theta_0_a))*(1-cos(delta_a))*cos(2*(alpha_a(l)-alpha_0_a))*sin(2*(alpha_a(l)-alpha_0_a))+sin(2*(theta_a(j)-theta_0_a))*(cos(delta_a)*cos(2*(alpha_a(l)-alpha_0_a))^2+sin(2*(alpha_a(l)-alpha_0_a))^2);
                        A(4) = cos(2*(theta_a(j)-theta_0_a))*sin(delta_a)*sin(2*(alpha_a(l)-alpha_0_a))-sin(2*(theta_a(j)-theta_0_a))*sin(delta_a)*cos(2*(alpha_a(l)-alpha_0_a));
                        
                        M = diag([1 1 -1 -1]);
                        
                        F(ind) = b(i)*A*M*G-I(ind);
                    end
                    
                    if (nargout > 1)
                        
                        if ~isnan(I(ind))
                            dG_dalpha0g(1) = 0;
                            dG_dalpha0g(2) = 2*(1-cos(delta_g))*sin(4*(alpha_g(k)-alpha_0_g)); % OK
                            dG_dalpha0g(3) = -2*(1-cos(delta_g))*cos(4*(alpha_g(k)-alpha_0_g)); % FIXED: the first factor was wrong
                            dG_dalpha0g(4) = 2*sin(delta_g)*cos(2*(alpha_g(k)-alpha_0_g));
                            
                            dG_dk1g(1) = 0;
                            dG_dk1g(2) = -sin(2*(alpha_g(k)-alpha_0_g))^2*sin(delta_g)/lambda(i); % FIXED: a "minus" sign was missing
                            dG_dk1g(3) = cos(2*(alpha_g(k)-alpha_0_g))*sin(2*(alpha_g(k)-alpha_0_g))*sin(delta_g)/lambda(i); % OK
                            dG_dk1g(4) = -sin(2*(alpha_g(k)-alpha_0_g))*cos(delta_g)/lambda(i);
                            
                            dG_dk2g = dG_dk1g/(lambda(i)^2); % OK
                            
                            dA_dalpha0a(1) = 0;
                            dA_dalpha0a(2) = 2*(cos(delta_a)-1)*(sin(2*(theta_a(j)-theta_0_a))*cos(4*(alpha_a(l)-alpha_0_a))-cos(2*(theta_a(j)-theta_0_a))*sin(4*(alpha_a(l)-alpha_0_a))); % OK
                            dA_dalpha0a(3) = 2*(cos(delta_a)-1)*(sin(2*(theta_a(j)-theta_0_a))*sin(4*(alpha_a(l)-alpha_0_a))+cos(2*(theta_a(j)-theta_0_a))*cos(4*(alpha_a(l)-alpha_0_a))); % OK
                            dA_dalpha0a(4) = -2*sin(delta_a)*(sin(2*(theta_a(j)-theta_0_a))*sin(2*(alpha_a(l)-alpha_0_a))+cos(2*(theta_a(j)-theta_0_a))*cos(2*(alpha_a(l)-alpha_0_a)));
                            
                            dA_dk1a(1) = 0;
                            dA_dk1a(2) = (-cos(2*(theta_a(j)-theta_0_a))*sin(2*(alpha_a(l)-alpha_0_a))^2+sin(2*(theta_a(j)-theta_0_a))*cos(2*(alpha_a(l)-alpha_0_a))*sin(2*(alpha_a(l)-alpha_0_a)))*sin(delta_a)/lambda(i); % OK
                            dA_dk1a(3) = (-sin(2*(theta_a(j)-theta_0_a))*cos(2*(alpha_a(l)-alpha_0_a))^2+cos(2*(theta_a(j)-theta_0_a))*cos(2*(alpha_a(l)-alpha_0_a))*sin(2*(alpha_a(l)-alpha_0_a)))*sin(delta_a)/lambda(i); % OK
                            dA_dk1a(4) = (cos(2*(theta_a(j)-theta_0_a))*sin(2*(alpha_a(l)-alpha_0_a))-sin(2*(theta_a(j)-theta_0_a))*cos(2*(alpha_a(l)-alpha_0_a)))*cos(delta_a)/lambda(i);
                            
                            dA_dk2a = dA_dk1a/(lambda(i)^2);
                            
                            dA_dtheta0a(1) = 0;
                            dA_dtheta0a(2) = 2*sin(2*(theta_a(j)-theta_0_a))*(cos(2*(alpha_a(l)-alpha_0_a))^2+cos(delta_a)*sin(2*(alpha_a(l)-alpha_0_a))^2)-2*cos(2*(theta_a(j)-theta_0_a))*(1-cos(delta_a))*cos(2*(alpha_a(l)-alpha_0_a))*sin(2*(alpha_a(l)-alpha_0_a)); % OK
                            dA_dtheta0a(3) = 2*sin(2*(theta_a(j)-theta_0_a))*(1-cos(delta_a))*cos(2*(alpha_a(l)-alpha_0_a))*sin(2*(alpha_a(l)-alpha_0_a))-2*cos(2*(theta_a(j)-theta_0_a))*(cos(delta_a)*cos(2*(alpha_a(l)-alpha_0_a))^2+sin(2*(alpha_a(l)-alpha_0_a))^2); % OK
                            dA_dtheta0a(4) = 2*sin(2*(theta_a(j)-theta_0_a))*sin(delta_a)*sin(2*(alpha_a(l)-alpha_0_a))+2*cos(2*(theta_a(j)-theta_0_a))*sin(delta_a)*cos(2*(alpha_a(l)-alpha_0_a));
                            
                            J(ind,1) = b(i)*A*M*dG_dalpha0g;
                            J(ind,2) = b(i)*A*M*dG_dk1g;
                            J(ind,3) = b(i)*A*M*dG_dk2g;
                            J(ind,4) = b(i)*dA_dalpha0a*M*G;
                            J(ind,5) = b(i)*dA_dk1a*M*G;
                            J(ind,6) = b(i)*dA_dk2a*M*G;
                            J(ind,7) = b(i)*dA_dtheta0a*M*G;
                            J(ind,8:end) = 0;
                            J(ind,7+i) = A*M*G;
                        end
                    end
                end
            end
        end
    end
        
end

