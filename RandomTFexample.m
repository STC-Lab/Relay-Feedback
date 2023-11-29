clear all, close all;
clc;
%% Another example with transfer function
Ts = .1;
s = tf('s');
G_ct = ((s+2.56)*(s-1.56))/((s+1)*(s+2)*(s+3));
G_dt = c2d(G_ct, Ts, 'tustin');
[num, den] = tfdata(G_ct, 'V');
% Do a partial fraction decomposition
[r, p, q] = residue(num, den);
i_neg = find(r<0); %find negative denominator
Gpos = 0;
%split into postive tranfser function and negative one
for i = 1:length(r)
    if i == i_neg
        Gneg = -tf(r(i),[1 -p(i)]);
    else
        Gpos =Gpos + tf(r(i),[1 -p(i)]);
    end
end
Gcheck = minreal(Gpos-Gneg);
%Check passivity of both tranfser functions
P_G = isPassive(G_dt);
P_Gneg = isPassive(Gneg);
P_Gpos = isPassive(Gpos);
%Check if the equality H = Hpos-Hneg is correct
[num_check,den_check] = tfdata(Gcheck, 'V');
B_num = isequaltol(num_check, num);
B_den = isequaltol(den_check, den);
fprintf(" Gneg is a passive transfer function: %s. \n Gpos is a passive transfer function: %s. \n Gcheck is equal to G_ct: %s %s \n", mat2str(P_Gneg),mat2str(P_Gpos), mat2str(B_num), mat2str(B_den))
%% Discretize newly found transferfunctions
Gd_neg = c2d(Gneg, Ts, 'tustin');
Gd_pos = c2d(Gpos, Ts, 'tustin');
SS_neg = ss(Gd_neg);
[Aneg,Bneg,Cneg,Dneg] = ssdata(SS_neg);
SS_pos = ss(Gd_pos);
[Apos,Bpos,Cpos,Dpos] = ssdata(SS_pos);
SS= ss(G_dt);
SS_ct=ss(G_ct);
[Act, Bct, Cct, Dct]  = ssdata(SS_ct);
[A,B,C,D] = ssdata(SS);
% Create Hankel matrices
[h, t] = impulse(G_dt);
H = hankel(h(2:end));
% Hneg = zeros(size(H));
% Hpos = zeros(size(H));
[hneg, tneg] = impulse(Gd_neg, t(end));
[hpos, tpos] = impulse(Gd_pos, t(end));

Hneg(1:length(hneg)-1,1:length(hneg)-1) = hankel(hneg(2:end));
Hpos(1:length(hpos)-1,1:length(hpos)-1) = hankel(hpos(2:end));

Hcheck = isequaltol(Hpos-Hneg, H);
%% Simulate system with relay feedback
k_sim = 2000000;
u = zeros(1,k_sim+1);
xpos = zeros(2,k_sim+1);
xpos(:,1) = randn(2, 1);
xneg = zeros(1,k_sim+1);
xneg(:,1) = randn(1, 1);
x = zeros(3,k_sim+1);
x(:,1) = randn(3, 1);
y = zeros(1,k_sim+1);
for i = 1:k_sim
    xpos(:,i+1) = Apos*xpos(:,i)+Bpos*u(i);
    xneg(:,i+1) = Aneg*xneg(:,i)+Bneg*u(i);
    x(:,i+1) = A*x(:,i)+B*u(i);
    y(i+1) = Cpos*xpos(:,i) - Cneg*xneg(:,i);
    u(i+1) = sign(y(i+1));
end
t_sim = (0:k_sim)*Ts;
figure;
plot(t_sim(1:100), u(1:100))
hold on
plot(t_sim(1:100), y(1:100))
xlabel('simulation time [s]')
ylabel('Relay output [-]')
hold off
figure
plot(xpos(1,:), xpos(2,:))
figure
plot3(x(1,:), x(2,:),x(3,:))
grid on
hold on
xyplane = polyshape([-60 -60 60 60],[60 -60 -60 60]);
pc_map = x(:,find(abs(x(3,:)) <3));
plot(xyplane)
scatter3(pc_map(1,:),pc_map(2,:),pc_map(3,:),'rx')
hold off
figure
scatter(pc_map(1,:),pc_map(2,:))
xlabel(x)
grid on
%% Poincaré mapping
x_eq = [0.6; -0.44; 0.32];
t_eq = 1.4;
v = Act*x_eq - Bct;
W = (eye(size(Act))-(v*Cct)/(Cct*v))*expm(Act*t_eq);
eig_W = eig(W);
if abs(eig_W) < 1
    fprintf("limitcycle is locally stable \n")
end
%simulate a trajectory
delta_x = zeros(3,k_sim/100000+1);
delta_x(:,1) = x_eq;
for i = 1:(k_sim/100000)
    delta_x(:,i+1) = W*delta_x(:,i);
end
% figure;
% plot(delta_x(1,:))
% hold on
% plot(delta_x(2,:))
% plot(delta_x(3,:))
% hold off

%% Checking the difference between various methods to determine the
%frequency response
% [Gd_pos_den, Gd_pos_num] = tfdata(Gd_pos,'V')
% FR_pos = ((Gd_pos_den(1)*j*k).^2+(Gd_pos_den(2)*j*k)+Gd_pos_den(3))./((Gd_pos_num(1)*j*k).^2+(Gd_pos_num(2)*j*k)+Gd_pos_num(3))% FR_2 =squeeze(FR_2)'
% for i = 1:length(k)
%     FR_1(i) = evalfr(Gd_pos,j*k(i));
% end
% FR_2 = freqresp(Gd_pos,j*k)
% FR_2 = flip(squeeze(FR_2))
% figure
% semilogx(k, abs(FR_pos), k,abs(FR_2) ,k,abs(FR_1))
% 
% isequaltol(FR_pos,FR_2, 0.1)
% isequaltol(FR_pos,FR_2, 0.1)
% isequaltol(FR_1,FR_2, 0.1)
% Difference was neglible, choose the easiest and foolproof method:
% freqresp
%% Limitcycle check
N = 500;
T = linspace(1,N,N);
ts = T(2)-T(1);
k = (1/ts)/N*(-N/2:N/2-1);
u_guess =[ones(1,N/2) -ones(1,N/2)]; %linspace(-1,1,N);
[x0, T] = zerofinding(Gpos,Gneg,T, N, k, u_guess);
legend
figure
plot(T, x0)
%% Functions
function [x0, T] = zerofinding(Gpos, Gneg,T, N, k, u_guess)
    epsilon = 0.0001;
    M = 1e6;
    alpha = 0.1;
    iters = 0;
    
    FRF_pos= flip(squeeze(freqresp((1+alpha*Gpos)^-1,2*pi*k)))';
    FRF_neg = flip(squeeze(freqresp(Gneg,2*pi*k)))';
    
    count = 0;
    x0 = u_guess;

    while true
        count = count +1;
        x1 = x0;
        x_star = lsim(Gneg,x0,T);
        x_star = x_star';
        x0 = DRsplitting_RF(FRF_pos, Gpos, x_star, x0);
        i = 1;
        for k = 1:length(x0)-1
            if x0(k)*x0(k+1) <= 0
                i = k;
                break
            end
        end
        %x0 = [x0(i:end) x0(1:i-1)];
        %check_out = max(abs(x0-x1))/max(abs(x0))
        if max(abs(x0-x1))/max(abs(x0)) < epsilon
            break
        elseif iters ~= 0 && count >= iters
            break
        elseif max(abs(x0-x1)) > M
            error("Unstable!!")
        elseif count == 1000
            break
        end
            fprintf("Outer iterations: %d\n", count)
            figure(6)

            plot(0:1:N-1,x0)
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


function i0 = DRsplitting_RF(FRF_pos, Gpos, x_star, x0)
    epsilon = 0.0001; M = 1e6; alpha = 0.1;
    count = 0;
    i0 = x0;
    while true
        count = count + 1;
        i1 = i0;
        y_fft = lsim((1+alpha*Gpos)^-1,i0+alpha*x_star,0:1:length(x0)-1);%%compute_resolvent(i0+alpha*x_star, FRF_pos);
        y2_fft =compute_resolvent(i0+alpha*x_star, FRF_pos); %compute_resolvent(x0,FRF_neg);
        % figure
        % plot(y_fft)
        % hold on
        % plot(y2_fft)
        % drawnow
        % legend('lsim', 'compute resolvent')

        %x_half =y2_fft;%
        x_half = y_fft';
        z_half = 2*x_half - i0;
        
        x = zeros(size(z_half));
        for i = 1:length(z_half)
            x(i) = resolvent_normalcone(z_half(i));
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
function  y = resolvent_normalcone(x)
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
function x = prox_l(v, l, u)
    lambda = 0.5;
    epsilon = 0.1;
    l = l+0.01;
    u = u -0.01;
    	x = v;
    while u-l > epsilon
        g = normalcone(x) + (1\lambda)*(x-v);
        a = x - lambda*g;
        b = x;
        if g <0
            a = x;
            b = x-lambda*g;
        end
        l = max(l ,a);
        u = min(u, b);
        x = (l+u)/2;
    end
end