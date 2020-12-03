% Theory ==================================================================
% initialise parameters
n1 = 1.0;
n2 = [1.0, 1.2, 1.4, 1.5, 1.7, 2.0, 3.0];
n3 = 1.5;
lam = 1.5e-6; %1500nm light
h = lam./ n2; %equal optical thicknesses for all different n2
theta2 = 0.0; %normal incidence

% evaluate for each n2
x = linspace(0, 1, 1000);
for i = 1:length(n2)
    H = linspace(0, h(i), length(x));
    y = [];
    for j = 1:length(x)
        y = [y ; R(n1, n2(i) , n3, lam , theta2, H(j))];
    end
    plot(x, y)
    hold on
end
chr = 'Reflectivity of a Dielectric Film of n_2';
chr = [chr newline 'as a Function of Optical Thickness'];
title(chr)
xlabel('Optical Thickness, λ_{0} [n_{2}h]')
ylabel('R')
hold off

% FUNCTIONS %==============================================================
% fresnel special case r and t coeff calculators
function value = rij(ni, nj)
    value = (ni - nj)/(ni + nj);
end

function value = tij(ni, nj)
    value = (2 * ni)/(ni + nj);
end

% reflection coeff of monoslab, r
function value = r(n1, n2, n3, lam0, theta2, h)
    i = sqrt(-1);
    beta = (2 * pi / lam0) * (n2 * h) * cos(theta2);
    r12 = rij(n1, n2);
    r23 = rij(n2, n3);
    numerator = (r12 + (r23 * exp(2 * i * beta)));
    denominator = (1 + (r12 * r23 * exp(2 * i * beta)));
    value = numerator / denominator;
end

% reflectivity of monoslab, R
function value = R(n1, n2, n3, lam0, theta2, h)
    beta = (2 * pi / lam0) * (n2 * h) * cos(theta2);
    r12 = rij(n1, n2);
    r23 = rij(n2, n3);
    numerator = (r12^2) + (r23^2) + (2 * r12 * r23 * cos(2 * beta));
    denominator = 1 + ((r12^2)*(r23^2)) + (2  * r12 * r23 * cos(2 * beta));
    value = numerator / denominator;
end