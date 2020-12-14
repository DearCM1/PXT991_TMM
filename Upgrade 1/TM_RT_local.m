% This code was last updated to check the results in comparision to the
% those obtained for a thin dielectric film in Born, Optics, page 64.

% Last Updated; 2020/12/02 20:18; Calum
% added generalisations to AB slab multipications for N=0.5,1.0,1.5,2.0...
% ======================================================== PARAMETERS === %
% first and last p_values in the TE polarisation
theta_1 = 0.0;
theta_2 = 0.0; %normal incidence
n_f = 1.00;
n_l = 1.50;
p_f = p(n_f, theta_1);
p_l = p(n_l, theta_1);

% arguments of the stratisfied media
n_1 = 1.45;
n_2 = 2.30;
d_1 = 120.7e-9;
d_2 = 76.1e-9;
p_1 = p(n_1, theta_1);
p_2 = p(n_2, theta_1);
N = 12;  %N periods of stratified media

% value of the PBG centre frequency
c = 2.99792458e+8;
D = d_1 + d_2;
n_avg = (d_1*n_1+d_2*n_2)/D;
w_0 = pi*c/(n_avg*D);

% range of frequencies
w_min = 0.5*w_0;
w_max = 1.5*w_0;

% ========================================================== TMM LOOP === %
%for loop to find R,T as a function of frequency
w_in = w_min:0.001*w_0:w_max;
T_out = zeros([1,length(w_in)]);
R_out = zeros([1,length(w_in)]);
for i = 1:length(w_in)
    w = w_in(i);
    M_A = TM(n_1, d_1, w, p_1);
    M_B = TM(n_2, d_2, w, p_2);
    %initialise and multiply the matrix AB layers N times
    M_N = [1 1;1 1];
    for j = 1:N*2
        if j == 1
            M_N = M_N .* M_B;
        elseif (-1)^j == -1
            M_N = M_N * M_B;
        else
            M_N = M_N * M_A;
        end
    end
    %reflection and trasmission coefficients
    a = (M_N(1,1) + M_N(1,2)*p_f - M_N(2,1)/p_l - M_N(2,2)*p_f/p_l);
    b = (M_N(1,2)*p_f + M_N(2,1)/p_l - M_N(2,2)*p_f/p_l - M_N(1,1));
    r = a / b;
    
    c = 2*((M_N(1,1)*M_N(2,2) - M_N(1,2)*M_N(2,1))*p_f/p_l);
    d = (M_N(1,1) + M_N(2,2)*p_f/p_l - M_N(1,2)*p_f - M_N(2,1)/p_l);
    t = c / d;
    %reflectivity and trasmittivity
    R_out(i) = abs(r).^2;
    T_out(i) = abs(t).^2*p_l/p_f;
end

% ========================================================== PLOTTING === %
% plot of R and log(T) vs w/w_0
rel_w = w_in/w_0;
figure(3)
plot(rel_w,R_out)
xlabel('w')
ylabel('R')
figure(2)
semilogy(rel_w,T_out)
xlabel('w/w_0')
ylabel('T')

% ============================================================ EXPORT === %
fileID = fopen('R_i.txt','w');
fprintf(fileID,'%6s %12s\n','R','w/w0');
fprintf(fileID,'%6.5f %12.5f\n',[R_out;rel_w]);
fclose(fileID);
fileID = fopen('T_i.txt','w');
fprintf(fileID,'%6s %12s\n','T','w/w0');
fprintf(fileID,'%6.5f %12.5f\n',[T_out;rel_w]);
fclose(fileID);

% ========================================================= FUNCTIONS === %
% p_value
function p_val = p(n, theta)
    p_val = n * cos(theta);
end

% transfer matrix
function M = TM(n,d,w,p)
    c = 2.99792458e+8;
    k = n*w/c;
    M = [cos(k*d) 1i*1/p*sin(k*d); 1i*p*sin(k*d) cos(k*d)];
end
