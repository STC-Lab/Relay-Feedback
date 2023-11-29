clear all,  close all;
clc;
%% Global variables
global alpha; %resolvent operator
global epsilon; % stopping treshold
global M; %Blowup treshold
alpha = 0.05;
epsilon = 0.01;
M = 1e6;
%% Choose static non linearity
% 1 = Relay, 2 = Saturation, 3 = Atanh
global Nonlinear;
Nonlinear = 1;
%% 
% Define Time
t_eq = 1;
N = 50000;
T_end = 30*t_eq;
t = linspace(0,T_end,N);
 
ts = t(2)-t(1);
%t = 0:ts:T_end;
%%
% Ts = 1;
% h = [1 1];
% H = hankel(h);
% a = ones(1,length(h)-1);
% A = diag(a,1);
% B = h';
% C = [1 zeros(1,length(h)-1)];
% D = 0;
% sys = ss(A,B,C,D,ts);
% Gh_d = tf(sys);
% Gh_c = d2c(Gh_d, 'Tustin');
% [num_h, den_h] = tfdata(Gh_c, 'V');
% num_h(abs(num_h) < 10e-2) = 0;
% num_h = round(num_h, 2);
% den_h(abs(den_h) < 10e-2) = 0;
% den_h = round(den_h, 2);
% Gh_c = tf(num_h, den_h);
%% Define Transfer function and do the splitting
s = tf('s');
G = -1*((s+2.56)*(s-1.56))/((s+1)*(s+2)*(s+3));
%G = (-s+1)/((s+1)*(s+2));
%G = minreal(Gh_c);
%G = -(s+7)/(s+5)^4;
%G = (2*s^2-4*s)/(s^2+4*s+4)
%Gd2 = c2d(G,.1, 'Tustin')
% Do a partial fraction decomposition
[num, den] = tfdata(G, 'v');
[r, p, q] = residue(num, den);
r = round(r);
p = round(p);
q = round(q);
i_neg = find(r<0); %find negative denominator
%find number of repeated poles

Gpos = 0;
Gneg = 0;
G_dummy = 0;
repeat_count = 0;
previous = 0;
%split into postive tranfser function and negative one
for i = 1:length(r)
    if i > 1     
        previous = isequal(p(i),p(i-1));
        if previous == 1;
            repeat_count = repeat_count +1
        else
            repeat_count = 0;
        end
    end
    if i > 1 && previous == 1
        if r(i) < 0
            G_dummy = r(i)/((s-p(i))^(repeat_count+1));
            Gneg =   Gneg - G_dummy;
        else
            G_dummy = r(i)/((s-p(i))^(repeat_count+1));
            Gpos = Gpos  + G_dummy;
        end
    else
        if r(i) < 0
            G_dummy = r(i)/((s-p(i)));
            Gneg = Gneg + -G_dummy;
        else
            G_dummy = r(i)/((s-p(i)));
            Gpos = Gpos + G_dummy;
        end
    end
end
if isempty(q) == 0
    if q > 0
            Gpos = Gpos +q;
    else
            Gneg = Gneg -q;
    end
end
Gneg = minreal(Gneg, 10e-2);
Gpos = minreal(Gpos, 10e-2);
Gcheck = minreal(Gpos-Gneg, 10e-2);
[num_Gcheck, den_Gcheck] = tfdata(Gcheck, 'V');
num_Gcheck(abs(num_Gcheck) < 10e-2) = 0;
den_Gcheck(abs(den_Gcheck) < 10e-2) = 0;
Gcheck = tf(num_Gcheck, den_Gcheck);
%Check passivity of both tranfser functions
P_Gneg = isPassive(Gneg);
P_Gpos = isPassive(Gpos);
%Check if the equality H = Hpos-Hneg is correct)
[num_check,den_check] = tfdata(Gcheck, 'V');
B_num = isequaltol(num_check, num);
B_den = isequaltol(den_check, den);
fprintf(" Gneg is a passive transfer function: %s. \n Gpos is a passive transfer function: %s. \n Gcheck is equal to G_ct: %s %s \n", mat2str(P_Gneg),mat2str(P_Gpos), mat2str(B_num), mat2str(B_den))
%% Make everyting discrete time
Gd = c2d(G,ts, 'Tustin');
Gd_neg = c2d(Gneg, ts, 'tustin');
Gd_pos = c2d(Gpos, ts, 'tustin');
P_G = isPassive(Gd);
P_Gneg = isPassive(Gd_neg);
P_Gpos = isPassive(Gd_pos);
%% Make state space
SS_neg = ss(Gd_neg);
[Aneg,Bneg,Cneg,Dneg] = ssdata(SS_neg);
SS_pos = ss(Gd_pos);
[Apos,Bpos,Cpos,Dpos] = ssdata(SS_pos);
SS= ss(Gd);
[A,B,C,D] = ssdata(SS);
SS_ct=ss(G);
[Act, Bct, Cct, Dct]  = ssdata(SS_ct);
% Hankel matrices
[hcheck, t_hankel] = impulse(Gd);
Hcheck = hankel(hcheck(2:end));
% Hneg = zeros(size(H));
% Hpos = zeros(size(H));
[hneg, tneg] = impulse(Gd_neg, t_hankel(end));
[hpos, tpos] = impulse(Gd_pos, t_hankel(end));

Hneg(1:length(hneg)-1,1:length(hneg)-1) = hankel(hneg(2:end));
Hpos(1:length(hpos)-1,1:length(hpos)-1) = hankel(hpos(2:end));

Hcheck = isequaltol(Hpos-Hneg, Hcheck);
%% Simulate system
k_sim = length(t);
u = zeros(size(t));
xpos = zeros(length(Bpos),k_sim);
xpos(:,1) = 5*randn(length(Bpos),1);
xneg = zeros(length(Bneg), k_sim);
xneg(:,1) = 5*randn(length(Bneg),1);
x = zeros(length(B),k_sim);
x(:,1) = randn(size(B));
y = zeros(size(t));
for i = 1:k_sim-1
    xpos(:,i+1) = Apos*xpos(:,i)+Bpos*u(i);
    xneg(:,i+1) = Aneg*xneg(:,i)+Bneg*u(i);
    %x(:,i+1) = A*x(:,i)+B*u(i);
    %y(i) = C*x(:,i);
    y(i) = Cpos*xpos(:,i) - Cneg*xneg(:,i);
    if Nonlinear == 1
        u(i+1) = -sign(y(i));
    elseif Nonlinear == 2
        u(i+1) = -saturation(y(i));
    elseif Nonlinear == 3
        u(i+1) = - atanh(y(i));
    end
end
figure
plot(t,u)
hold on
zci = @(v) find(diff(sign(v)));
i_check = zci(u);
scatter(t(i_check),u(i_check), 'x');
for i = 1:length(i_check)-1
    period_check(i) = t(i_check(i+1))-t(i_check(i));
end
% figure
% plot(t,y)
%% Poincaré mapping
%Check stability of limitcyle
% x_eq = [0.6; -0.44; 0.32];
% v = Act*x_eq - Bct;
% W = (eye(size(Act))-(v*Cct)/(Cct*v))*expm(Act*t_eq);
% eig_W = eig(W);
% if abs(eig_W) < 1
%     fprintf("limitcycle is locally stable \n")
% end
% % Make poincaré map
% x_pm = zeros(length(B),k_sim);
% x_pm(:,1) = x_eq;
% g = -expm(Act*t_eq)*x + expm(Act*t_eq)*Act^-1*Bct;
% figure
% plot3(x(1,:),x(2,:),x(3,:))
% grid on
% hold on
% scatter3(g(1,:),g(2,:),g(3,:))
%% Zero finding algorithm

u_guess = square((0.25*pi)*t);

%u_guess = linspace(-1,1,length(t));
Fs = 1/ts;
k = Fs/N*(-N/2:N/2-1);
[x0, T] = zerofinding(Gpos, Gneg, Gd_pos, Gd_neg, t, N, k, u_guess);
figure
plot(T, x0)
hold on
i_zc = zci(x0);
scatter(T(i_zc),x0(i_zc), 'x')
for i = 1:length(i_zc)-1
    period(i) = T(i_zc(i+1))-T(i_zc(i));
end
[pks, locs] = findpeaks(x0,'MinPeakProminence',0.2);
%save("Hankel[1 1 0]")

%% Functions


function [x0, T] = zerofinding(Gpos, Gneg, Gd_pos, Gd_neg, T, N, k, u_guess)
    global epsilon;
    global M;
    global alpha;
    iters = 0;
    
    FRF_pos= flip(squeeze(freqresp((1+alpha*Gd_pos)^-1,2*pi*k)))';
    FRF_neg = flip(squeeze(freqresp(Gd_neg,2*pi*k)))';
    
    count = 0;
    x0 = u_guess;

    while true
        count = count +1;
        x1 = x0;
        %x2_star = lsim(Gneg, x0, T);
        %x_star = x2_star';
        x_star = compute_resolvent(x0, FRF_neg);
        x0 = DRsplitting_RF(Gpos, Gd_pos, FRF_pos, x_star, x0, T);
        % phase shift the time axis
        i = 1;
        for k = 1:length(x0)-1
            if x0(k)*x0(k+1) <= 0
                i = k;
                break
            end
        end
        x0 = [x0(i:end) x0(1:i-1)];


        %check_out = max(abs(x0-x1))/max(abs(x0))
        if max(abs(x0-x1))/max(abs(x0)) < epsilon
            break
        elseif iters ~= 0 && count >= iters
            break
        elseif max(abs(x0-x1)) > M
            error("Unstable!!")
        elseif count == 3000
            break
        end
            fprintf("Outer iterations: %d\n", count)
            figure(6)
            hold off
            plot(T,x0)
            
            drawnow
            %pause(2)
        if false
            i = 1;
            for k = length(x0)-1:-1:1
                if x0(k)*x0(k+1) <= 0
                    i = k;
                    break
                end
            end
            offset = length(x0)-i*1/N;
            period = 1 - offset/2;
            T = linspace(0, N, period);
        end
    end
end


function i0 = DRsplitting_RF(Gpos, Gd_pos, FRF_pos, x_star, x0, T)
    global epsilon; global M; global alpha;
    global Nonlinear;
    count = 0;
    i0 = x0;
    while true
        count = count + 1;
        i1 = i0;
        y_fft = lsim((1+alpha*Gpos)^-1,i0+alpha*x_star,T);%%
        
        y2_fft =compute_resolvent(i0+alpha*x_star, FRF_pos); %compute_resolvent(x0,FRF_neg);


        x_half =y2_fft;
        %x_half = y_fft';
        z_half = 2*x_half - i0;
        
        x = zeros(size(z_half));
        for i = 1:length(z_half)
            if Nonlinear == 1
                %x(i) = saturation(z_half(i));
                x(i) = prox_l(z_half(i), -1, 1);
            elseif Nonlinear == 2
                x(i) = resolvent_saturation(z_half(i), alpha);
            elseif Nonlinear == 3
            x(i) = prox_gN(z_half(i), -1, 1);
            end
        end


        i0 = i0 + x - x_half;
        % figure(7)
        %     if mod(count,10) == 0
        %     plot(1:1:length(x0),i0), hold on
        %     drawnow
        %     end
        check = max(abs(i0-i1))/max(abs(i0));
        if max(abs(i0-i1))/max(abs(i0))< epsilon
            break
        elseif max(abs(i0-i1))> M
            error("Unstable!!")
        end
        if mod(count,100) == 0
            fprintf("Inner iterations : %d\n", count)
        end
        %pause(2)
        
    end
end


function y_t = compute_resolvent(u, sys_fr)
    uf = fftshift(fft(u));
    yf = sys_fr.*uf;
    y_if = ifft(ifftshift(yf));
    y_t = real(y_if);
end

function y = normalcone(u)
    for i = 1:length(u)
        if abs(u(i)) < 1
            y(i) = 0;
        else
            y(i) = NaN;
        end
    end
end
function  y = saturation(x)
    for i = 1:length(x)
        if abs(x) < 1
            y= x;
        elseif x>= 1
            y = 1;
        elseif x <= -1
            y =-1;
        end
    end
end

function y = resolvent_saturation(x, alpha)
    if abs(x)< 1+alpha
        y = 1/(1+alpha)*x;
    elseif x >= 1+alpha
        y = 1;
    elseif x <= -1-alpha
        y = -1;
    end
end

function y = resolvent_relu(x, alpha)
    if x < 0
        y = 0;
    else
        y = (1+alpha)^(-1)*x;
    end
end

function x = prox_l(v, l, u)
   global alpha;
   global epsilon;
   x = v;
   while u - l > epsilon
       g = normalcone(x) + (1/alpha)*(x - v);
       a = x - alpha*g;
       b = x;
       if g< 0
           a = x;
           b = x- alpha*g;
       end
       l = max(l, a);
       u = min(u, b);
       x = (l+u)/2;
   end
end

function x = prox_gN(v, l, u)
    global alpha;
    mu = 0.01;
    lambda = 0.5; %guard parameter

    x = v;
    if v< l
        x = l + mu;
    elseif v > u
        x = u - mu;
    end


   g = atanh(x)+(1/alpha)*(x-v);
   a = x - alpha*g;
   b = x;
   if g <0
       a = x;
       b = x - alpha*g;
   end
   
   l = max(l,a);
   u = min(u, b);
   x = (l+u)/2;
   while u-l > mu
       al = (u+l)/2 - lambda*(u-l)/2;
       au = (u+l)/2 + lambda*(u-l)/2;

       df = atanh(x) + (1/alpha)*(x-v);
       ddf = dtanhinv(x) + 1/alpha;
       update = x - df/ddf;
       
       if update < al 
           x = al;
       elseif x > au 
           x = au;
       end

       g = atanh(x) + (1/alpha)*(x-v);
       a = x - alpha*g;
       b = x;
       if g<0
           a = x;
           b = x-alpha*g;
       end
       l = max(l, a);
       u = min(u, b);
       x = (l+u)/2;
    end
end

function y = dtanhinv(x)
    if abs(x) < 1
        y = 1/(1-x^2);
    else
        y = NaN;
    end
end

