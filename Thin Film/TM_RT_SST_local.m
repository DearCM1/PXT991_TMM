% =====================================================TMM PARAMETERS === %
% first and last p_values in the TE polarisation
theta_1 = 0.0;
theta_2 = 0.0; %normal incidence
n_f = 1.00;
n_l = 1.50;
p_f = p(n_f, theta_1);
p_l = p(n_l, theta_1);

% arguments of the stratisfied media
n_val = 6; % what n_1 value to be used in calculations
n_1 = [1.0, 1.2, 1.4, 1.5, 1.7, 2.0, 3.0];
n_2 = 0.00; % not in use for thin film
d_1 = 0.00; % calculated later
d_2 = 0.00; % not in use for thin film
p_1 = p(n_1, theta_1);
p_2 = p(n_2, theta_1); % not in use for thin film
N = 0.5;  %N periods of stratified media = 1 layer

% initial calculations
lambda = 1.500e-6;
h = lambda./ n_1; % equal optical thicknesses for all different n_1
d_1 = h;
w_0 = w_calc(lambda);

% ========================================================== TMM LOOP === %
% for loop to find R,T as a function of film thickness
x = linspace(0, 1, 1000);
H = linspace(0, h(n_val), length(x));
T_out = zeros([1,length(x)]);
R_out = zeros([1,length(x)]);
for i = 1:length(H)
    M_A = TM(n_1(n_val), H(i), w_0, p_1(n_val));
    M_B = TM(n_2, d_2, w_0, p_2);
    %initialise and multiply the matrix AB layers N times
    M_N = [1 1;1 1];
    for j = 1:N*2
        if j == 1
            M_N = M_N .* M_A;
        elseif (-1)^j == -1
            M_N = M_N * M_B;
        else
            M_N = M_N * M_A;
        end
    end
    %reflection and trasmission coefficients
    a = (M_N(1,1) + M_N(1,2)*p_f - M_N(2,1)/p_l - M_N(2,2)*p_f/p_l);
    b = (M_N(1,2)*p_f + M_N(2,1)/p_l - M_N(2,2)*p_f/p_l - M_N(1,1));
    ri = a / b;
    
    c = 2*((M_N(1,1)*M_N(2,2) - M_N(1,2)*M_N(2,1))*p_f/p_l);
    d = (M_N(1,1) + M_N(2,2)*p_f/p_l - M_N(1,2)*p_f - M_N(2,1)/p_l);
    ti = c / d;
    %reflectivity and trasmittivity
    R_out(i) = abs(ri).^2;
    T_out(i) = abs(ti).^2*p_l/p_f;
end

% ============================================================ THEORY === %
% evaluate R theory function in same range
y = linspace(0,0,length(x));
for i = 1:length(x)
    y(i) = R(n_f, n_1(n_val) , n_l, lambda, theta_2, H(i));
end

% ============================================================ EXPORT === %
fileID = fopen('R_i.txt','w');
fprintf(fileID,'%6s %12s %12s\n','H','R_TMM','R_Theory');
fprintf(fileID,'%6.4f %12.4f %12.4f\n',[x;R_out;y]);
fclose(fileID);

% ========================================================== PLOTTING === %
figure(1)
plot(x, y)
hold on
plot(x, R_out)
hold off
xlim([0,1])
ylim([0,1])
xlabel('h/(lambda_0/n_1)')
ylabel('R')
legend('Theory','TMM')

figure(2)
plot(x, T_out)
xlim([0,1])
ylim([0,1])
xlabel('h/(lambda_0/n_1)')
ylabel('T')
legend('TMM')

% ========================================================= FUNCTIONS === %
% TMM =====================================================================
% p_value
function p_val = p(n, theta)
    p_val = n * cos(theta);
end

%transfer matrix
function M = TM(n,d,w,p)
    c = 3e+8;
    k = n*w/c;
    M = [cos(k*d) 1i*1/p*sin(k*d); 1i*p*sin(k*d) cos(k*d)];
end

% Theory ==================================================================
% angular frequency calculator
function ang_freq = w_calc(lambda)
    c = 2.99792458e+8;
    ang_freq = 2 * pi * c / lambda;
end

% fresnel special case r and t coeff calculators
function value = rij(ni, nj)
    value = (ni - nj)/(ni + nj);
end

%reflection coeff of monoslab, r
function value = r(n1, n2, n3, lam0, theta_2, h)
    i = sqrt(-1);
    beta = (2 * pi / lam0) * (n2 * h) * cos(theta_2);
    r12 = rij(n1, n2);
    r23 = rij(n2, n3);
    numerator = (r12 + (r23 * exp(2 * i * beta)));
    denominator = (1 + (r12 * r23 * exp(2 * i * beta)));
    value = numerator / denominator;
end

%reflectivity of monoslab, R
function value = R(n1, n2, n3, lam0, theta_2, h)
    beta = (2 * pi / lam0) * (n2 * h) * cos(theta_2);
    r12 = rij(n1, n2);
    r23 = rij(n2, n3);
    numerator = (r12^2) + (r23^2) + (2 * r12 * r23 * cos(2 * beta));
    denominator = 1 + ((r12^2)*(r23^2)) + (2  * r12 * r23 * cos(2 * beta));
    value = numerator / denominator;
end