close all, clear all;
clc;
%% Generate impulse response
Ts = 1;
T_end = 20;
p = 0.8;
t = 0:Ts:T_end;
h = zeros(1, length(t));
for i = 1:length(h)
    if t(i) < 2
        h(i) = t(i)-1;
    elseif t(i) >= 2 && t(i) < 2+p
        h(i) = 1-(t(i)-2)/p;
    else
        h(i) = 0;
    end
end
t_stop = find(t>2+p,1, 'first');
figure
plot(t,h)
%% Create state space

a = ones(1,length(h)-1);
A = diag(a,1);
B = h';
C = [1 zeros(1,length(h)-1)];
D = 0;
sys = ss(A,B,C,D,Ts)
G = tf(sys)
[Z, P, K] = zpkdata(sys, 'Vector');
% find zeroes outside of unit circle
Z_nmp = Z(abs(Z) > 1)
%% Creating Hankel matrix and splitting in Hp and Hn
H = hankel(h);
hneg = min(h, 0)
hpos = max(h, 0)
H1 = hankel(hpos)
H2 = hankel(hneg)
% Based on Eigendecompositon
[M, L] = eig(H);
[row, col] = find(L < 0);
Lneg = zeros(length(h), length(h));
Lneg(:,col) = L(:,col);
Hneg = M*sqrt(Lneg)*(M*sqrt(Lneg))';
Hpos = H + Hneg;
Hcheck = (Hpos-Hneg);
check = isequal(H, Hcheck);
fprintf('The condition Hcheck is equal to H is %s\n', mat2str(check))
% Alternating antidiagonals
h_alt1 = zeros(1,length(h))
h_alt1(1:2:end) = h(1:2:end)
H_alt1 = hankel(h_alt1)
h_alt2 = zeros(1,length(h))
h_alt2(2:2:end) = h(2:2:end)
H_alt2 = hankel(h_alt2)

%Regular matrix splitting
Hdiag = diag(diag(H))
Hnodiag = H - H1;
%% system with relay feedback response
k_sim = 100;
u = zeros(1,k_sim+1);
x = zeros(length(h),k_sim+1);
x(:,1) = randn(length(h), 1)
y = zeros(1,k_sim+1);
for i = 1:k_sim
    x(:,i+1) = A*x(:,i)+B*u(i);
    y(i+1) = C*x(:,i);
    u(i+1) = sign(y(i+1));
end
figure;
plot(0:k_sim, u)
hold on
plot(0:k_sim,y)

