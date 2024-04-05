close all, clear all;
clc;
%% choose transfer function
s = tf('s');
%G = 1/(s^3+5*s^2+2*s+1);
%G = (-s+1)/((s+1)*(s+2));
%G = (s+1)/((s+1)*(s+2));
G = -1*((s+2.56)*(s-1.56))/((s+1)*(s+2)*(s+3));
% minphase
%G = -1*((s+2.56)*(s+1.56))/((s+1)*(s+2)*(s+3));
%G = exp(-s)*s^2/(s^3+2*s^2+2*s+1)
%% Make Nyquist plot of transfer function and calculate the Beltrami-Klein convex hull to determine the SRG
w = logspace(-2, 3, 300);
[r, i] = nyquist(G, w);
r = squeeze(r);
i = squeeze(i);
nyq = r +j*abs(i);
nyq_p = BeltramiKlein_PW(nyq);
k = convhull([real(nyq_p) imag(nyq_p)]);
srg_k = nyq_p(k);
srg_k = srg_k(1:end-1);
srg = invBeltramiKlein_PW(srg_k)';
arc = arc_min(srg(1), srg(end), 50);
srg = [srg flip(arc)];
srg_full = [srg (srg').'];
%% Plots
figure('Name','SRG and Nyquist')
fill(real(srg_full), imag(srg_full), [0.83 0.83 0.83], 'LineStyle','none')
hold on
n = nyquistplot(G);
nyqopt = getoptions(n);
nyqopt.Title.String = 'Nyquist plot and SRG of system G(s)';
setoptions(n,nyqopt);
h = findall(gcf, 'Type', 'Line');
for i = 4:length(h)         % Loop through every patch object handle.
    h(i).LineWidth = 2;     % Set the new LineWidth value.
end
l = findall(gcf, 'Type', 'Patch');
for i = 1:length(l)         % Loop through every patch object handle.
    l(i).LineWidth = 3;     % Set the new LineWidth value.
end
grid off
rectangle('Position', [-1 -1 2 2], 'Curvature', [1 1], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle','--');


figure('Name', 'SRG inverse')
srgfull_inv = complex_inv(srg_full);
fill(real(srgfull_inv), imag(srgfull_inv), [0.83 0.83 0.83], 'LineStyle', 'none')
hold on
rectangle('Position', [-1 -1 2 2], 'Curvature', [1 1], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle','--');
title("SRG of G^{-1}")
xlim([-5 5]), ylim([-5 5])
plot([0 0], ylim, ':', 'Color', [0 0 0])
plot(xlim, [0 0], ':', 'Color', [0 0 0])
xlabel('Real Axis'), ylabel('Imaginary Axis')

figure("Name",'SRG inverse zoomed out')
fill(real(srgfull_inv), imag(srgfull_inv), [0.83 0.83 0.83], 'LineStyle', 'none')
hold on
title("SRG of G^{-1}")
plot([0 0], ylim, ':', 'Color', [0 0 0])
plot(xlim, [0 0], ':', 'Color', [0 0 0])
xlabel('Real Axis'), ylabel('Imaginary Axis')

%% SRG of static nonlinearity
h = 0.001;
x = -10:0.001:10;
rel = sign(x);
satv = sat(x);

SRG_rel = SRG_Nonlin(rel, h);
SRG_sat = SRG_Nonlin(satv, h);

%% Feedback interconnection (assume G is incrementally stable)
% invert srg of G 
srg_inv = complex_inv(srg_full);

% make Polygon of the srg's
SRGinv_pgon = polyshape(real(srg_inv), imag(srg_inv));
Rel_pgon = polyshape(real(SRG_rel), imag(SRG_rel));
Sat_pgon = polyshape(real(SRG_sat), imag(SRG_sat));


% Minkowski sum of G inverse plus rleay

GRel_pgon = minkowskiSum(SRGinv_pgon, Rel_pgon);
GSat_pgon = minkowskiSum(SRGinv_pgon, Sat_pgon);
figure('Name','Minkowski sum relay')
plot(GRel_pgon);

figure('Name','Minkowski sum saturation')
plot(GSat_pgon);

%% Approximate SRG of relay feedback system
[re_GRel, im_GRel] = boundary(GRel_pgon);


% Left boudnary of SRG sum
x_left = -5*ones(1,20000);
y_left = linspace (-10000,10000,20000);
z_left = x_left +j*y_left;
zinv_left = complex_inv(z_left);

% grid of uniformly distributed points inside SRG
%xmin = min(re_GRel);xmax = max(re_GRel); ymin = min(im_Grel); ymax = max(im_Grel);
xmin = -5; xmax = 50; ymin = -100; ymax = 100;
rect = polyshape([xmin xmin xmax xmax], [ymin ymax ymax ymin]);
x_grid = (xmax-xmin).*rand(10000,1)+xmin;
y_grid = (ymax-ymin).*rand(10000,1)+ymin;
in_shape = inpolygon(x_grid, y_grid, re_GRel, im_GRel);
x_grid = x_grid(in_shape); y_grid = y_grid(in_shape);
z_grid = x_grid+j*y_grid;
srg_RFS = complex_inv(z_grid);
%figure, scatter(x_grid, y_grid);
% sum_complex = real_sum+j*imag_sum;
% srg_RFS = complex_inv(sum_complex);
%% Plot SRG of RFS
figure('Name','SRG Relay Feedback system');
scatter(real(srg_RFS), imag(srg_RFS));
hold on
plot(zinv_left)
grid on
xlim([-5 5]), ylim([-5 5])
%% SRG of saturation feedback system
[re_GSat, im_GSat] = boundary(GSat_pgon);
complex_Gsat = re_GSat+j*im_GSat;
SRG_SFS = complex_inv(complex_Gsat);
SRG_SFSpgon = polyshape(real(SRG_SFS), imag(SRG_SFS));
figure('Name', 'SRG saturation feedback')
plot(SRG_SFSpgon, "LineStyle","none", "FaceColor",[0.6 0.6 0.6]);
%fill(real(SRG_SFS), imag(SRG_SFS), [0.83 0.83 0.83], 'LineStyle', 'none')
hold on
rectangle('Position', [-1 -1 2 2], 'Curvature', [1 1], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle','--');
plot([0 0], ylim, ':', 'Color', [0 0 0])
plot(xlim, [0 0], ':', 'Color', [0 0 0])
xlabel('Real Axis'), ylabel('Imaginary Axis')
%% Overapproximation of Saturation feedback system
half_cirlce = arc_min(j, -j, 350);
i_axis = linspace(-j,j, 350);
srg_SAT = [half_cirlce i_axis];

% Calculate minimum distance between saturation and Ginverse
distance = [];
for i = 1:length(srg_inv)
    for k = 1:length(srg_SAT)
        distance(i,k) = norm(srg_inv(i)-srg_SAT(k));
    end
end
d_min= min(distance(:));
[row, col] = find(distance == d_min);
points_G = srg_inv(row(1));
points_Sat = srg_SAT(col(1));
dmin_line = [points_G points_Sat];
dline_x = real(dmin_line); dline_y = imag(dmin_line);
% Plot SRG of G^-1 and SAT
figure('Name','G^{-1} and Saturation')
fill(real(srgfull_inv), imag(srgfull_inv), [0.83 0.83 0.83], 'LineStyle', 'none')
hold on
fill(real(srg_SAT), imag(srg_SAT), [0.6 0.6 0.6], 'LineStyle','none')
rectangle('Position', [-1 -1 2 2], 'Curvature', [1 1], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle','--');
xlim([-2 2]), ylim([-2 2])
plot(dline_x,dline_y, 'r--', 'LineWidth', 1.5)
% arrow = annotation('doublearrow')
% arrow.Parent = gca;
% arrow.Position = [dline_x(2), dline_y(2), dline_x(1), dline_y(1)];
plot([0 0], ylim, ':', 'Color', [0 0 0])
plot(xlim, [0 0], ':', 'Color', [0 0 0])
xlabel('Real Axis'), ylabel('Imaginary Axis')

%plot of SRG
figure('Name','SRG Saturation feedback')
rectangle('Position',[-1/d_min -1/d_min 2*1/d_min 2*1/d_min], 'Curvature',[1 1], 'EdgeColor','none', 'FaceColor',[0.83 0.83 0.83])
hold on
plot([0 0], ylim, ':', 'Color', [0 0 0])
plot(xlim, [0 0], ':', 'Color', [0 0 0])
xlabel('Real Axis'), ylabel('Imaginary Axis')

annotation('textarrow', [0.52 0.7750], [0.52 0.8150], 'String', '1/r_m')
%% Functions
function Y = BeltramiKlein_Op(A)
    n = length(A);
    Y = (I+A'.*A)^(-0.5).*(A'-j*I)*(A-j*I).*(I+A'.*A).^(-0.5);
end

function y = BeltramiKlein_PW(Z)
    y = zeros(size(Z));
    for i = 1:length(Z)
        z = Z(i);
        y(i) = ((conj(z)- j)*(z-j))/(1+conj(z)*z); 
    end
end

function Z = invBeltramiKlein_PW(f)
    Z = zeros(size(f));
    for i = 1:length(f)
        fz = f(i);
        Z(i) = (imag(fz)+j*sqrt(1-abs(fz)^2))/(real(fz)-1);
    end
end

% function c = ccw(x, y, z)
%     if angle(y - x) < angle(z - x)
%         c = 1;
%     elseif angle(y - x) > angle(z - x)
%         c = -1;
%     else
%         c = 0;
%     end
% end

function arcmin = arc_min(z1, z2,n) %function to create the smallest circular arc between two imaginary numbers
    x1 = real(z1);
    x2 = real(z2);
    y1 = imag(z1);
    y2 = imag(z2);
    
    if z1 == z2;
        arcmin = 0;
        return
    elseif x1-x2 ~= 0
        xc = (y1 - y2)*(y1 + y2)/2/(x1 - x2) + (x1 + x2)/2;
    else
        xc = (x1 - x2)*(x1 + x2)/2/(y1 - y2) + (y1 + y2)/2;
    end

    r = sqrt(x1-xc)^2+y1^2;

    theta1 = 0;
    theta2 = 0;
    theta1 = angle(z1 - xc + j*0);
    theta2 = angle(z2 - xc + j*0);

    arcmin = [];
    phi = linspace(theta1,theta2,n);
    for i = 1:length(phi)
        arcmin = [arcmin xc+r*exp(j*phi(i))];
    end
end
        
function z_inv = complex_inv(z) %% Möbius inverse, while keeping the sign of the angle the same, i.e. points inside the unit cirlce map to outside the unit circle and vice versa
    r = abs(z);
    a = angle(z);
    z_inv = 1./r.*exp(j.*a);
end
% function to overapproximate SRG of sector bouned nonlinearity 
function srg = SRG_Nonlin(fun, h)
    n = 100;% number of data points per line segment
    df = diff(fun)/h;
    df_max = max(df);
    df_min = min(df);
    if df_min >= 0 %% check if function is monotone
        if df_max <= 100 %check if function is Lipschitz, some arbitrary bound on the maximum slope, this is dependent on the sampling of your signal and your preference
            srg = arc_min(df_max*j, -df_max*1j, n); %right half of a circle centered at the origin with radius df_max
            srg = [srg linspace(-df_max*j, df_max*j, n)]; %complete the half circle with the imaginary axis
        else
            srg = [linspace(df_min-1000*j,df_min+1000*j,n) linspace(df_min+1000*j, (df_min+1000)+1000*j, n) linspace((df_min+1000)+1000*j, (df_min+1000)-1000*j, n) linspace((df_min+1000)-1000*j,df_min-1000*j, n)]; %section of the right half plane
        end
    elseif df_max <= 100 %check if function is Lipschitz, some arbitrary bound on the maximum slope, this is dependent on the sampling of your signal and your preference
        srg = [arc_min(df_max*j, -df_max*1j, n) -arc_min(df_max*j, -df_max*1j, n)]; %circle with centered at the origin with radisu df_max
    else
        fprintf('SRG cannot be approximated \n')
        return;
    end
end

function  y = sat(x)
    for i = 1:length(x)
        if abs(x(i)) < 1
            y(i)= x(i);
        elseif x(i)>= 1
            y(i) = 1;
        elseif x(i) <= -1
            y(i) = -1;
        end
    end
end

function y = somnonlin(x)
    for i = 1:length(x)
        if abs(x(i)) < 1
            y(i)= 50*x(i);
        elseif x(i)>= 1
            y(i) = 100*x(i);
        elseif x(i) <= -1
            y(i) = 100*x(i);
        end
    end
end
    
    