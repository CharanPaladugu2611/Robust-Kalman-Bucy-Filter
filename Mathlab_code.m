%% Robust Kalman-Bucy Filter
% Article author: Jemin George
% ENPM667 - Controls
% Students in team: Jason Chen and Venkata Charan Sai Palagdugu

%% Common Variables
A = [0, 0, -1, 0; ...
    0, 0, 0, -2; ...
    2.3, 1.3, -8.5, -1.62; ...
    1.78, 2.03, -1.02, -9.5];

C = [4.3, 0, -2.31,  0; ...
    0, -2.4 , 0, 1.54; ...
    3.64,  2.74, 0, 0; ...
    0, 2.42, 1.12, 4.34];

Am = [0, 0, -1, 0; ...
    0, 0, 0, -2; ...
    .3, .13, -4.24, .2; ...
    -.4, .53, .42, -5.43];

Cm = [.43, 0, -.231, 0; ...
    0, -2.4, 0, 1.54; ...
    .95, 1.7, 0, 0; ...
    0, 0, 1.3, .43];

% 6x6 matrix
F = [[Am, zeros(4,2)];
    [zeros(2,4), -eye(2,2)]];
    
% 4x6 matrix
H = [Cm,  [zeros(2,2); eye(2,2)]];
Gm= [eye(2,2); 
    zeros(2,4).'];
X01 = [1; -2; 5; 0]*ones(1,4);
X02 = [1; -2; 5; 0];
Xm0 = zeros(4,1);

Q = [0.25, 0; 0, 0.35];

% We have to modify the input values a little to make the Ricatti equation
% function accurately:
GQG = Gm*Q*Gm.';

% the time span used in the paper:
tv = linspace(0,40,200).'; 

%% Simulation of Inputs and Noise
close all
W1 = -5*diric(tv,8) - (tv/40).*exp(-tv/40);
W2 = 3*square(tv*2*pi/14,50) - (tv/40).*exp(-tv/40);
figure;
plot(tv,W1);
grid on;
figure;
plot(tv,W2);
grid on;
R = eye(4,4)*1E-4;
% we have to modify the input values a little to make the Ricatti equation
% work correctly:
[P1,K1,~] = icare(F.',H.',GQG,R,[],[],[]); % RICATTI EQUATION BABYYYYYY
K1 = K1.'
% K1_1 = P1*H.'*inv(R);
% K1 found in both ways are pretty much the same, but for the former
% option, we need to tranpose the matrix using .' then it is the exact same.

%% Matrix Inquality Derivation
close all
S = ones(6,6); % S is positive-definite matrix, but we do not know the values
setlmis([]); % begin setting up linear matrix inequalities
FKH = (F - K1*H); % crucial term, often reused!
E = [zeros(2,4); eye(4,4)]; % 6x4 matrix
X_ = lmivar(2,[6,6]); % 6x4 matrix
A_ = lmivar(2,[6,4]); % 6x4 matrix
% the 's' flag indicates symmetry, this is the only way that this LMI can
% be feasible even though it is not accurate
lmiterm([1 1 1 X_], FKH.', 1); % top left term (1,1), not accurate, but allows for feasibility
lmiterm([1 1 1 0], S); % top left   term (1,1)
lmiterm([1 1 1 X_], 1, FKH); % top left  term (1,1)
lmiterm([1 1 2 X_ ], 1, 1); % top right term (1,1)
lmiterm([1 1 2 -A_], -H.', 1); % top right term (1,1)
% lmiterm([1 2 1 -X_], 1, 1); % bottom left term (1,1)
% lmiterm([1 2 1 A_], -1, H); % bottom left term (1,1)
lmiterm([-1 2 2 0],0); % bottom right term (1,1), not necessary since
% MATLAB defaults to 0 anyways
lmis = getlmis();
[tmin,xfeas] = feasp(lmis);
aX_ = dec2mat(lmis,xfeas,X_) % 4x6 matrix
aA_ = dec2mat(lmis,xfeas,A_) % 4x6 matrix

%% Error Estimation 
close all
disp("==================== Complete! ============================")
% In this case, A = Am, so this cancels out into a zero matrix
dA = zeros(4);
dC = zeros(4);
psi = [dA; dC*(Am + eye(4,4) + dA)]; % 8x4 matrix
% New equation is adE = aF*aE + g*aB
% new state variables are aZmt, aZt, and aXt
% expanded noise is given by aB = -[dVt, dBt]
% We are looking to solve for matrix P a.k.a. covariances matrix

% 99% sure these values are correct, results in no errors, but I had to change a
% few things in the paper to make it work...
n = 4; 
m = 2;

% g is a 16x6 matrix
g = [-K1,                Gm;         % this row is a 6x6 matrix 
        zeros(n,m),  eye(n,n); % this row is a 4x6 matrix
        -K1,               Gm];        % this row is a 6x6 matrix
% aA_2 = [-16.0721, -0.5073, 0.1358, 8.0142; ...
%                 0.0098, 7.7899, 0.5412, 2.0673; ...
%                 -6.0882, 3.3429, 12.1134, 10.7137;...
%                 16.5996, -8.6149, -10.7137, 35.6611];

% 99% sure this is 16x6 matrix
aE = eye(16,6); % variable containing the new state variables

omg = [eye(n+m,n+m), zeros(n+m,n), eye(n+m,n+m)]; %this row is 6x16

% J should be a 6x6 matrix if the other matrices are true-to-size
J = norm(aA_*H*omg*aE)^(-3) * ...
    (-(aA_*H*omg*aE).'*aA_*H*omg*aE*eye(n+m) + ...
    aA_*H*omg*aE*(aA_*H*omg*aE).')*aA_*H;

% aF is a x matrix
col1 = [0*(FKH + J); zeros(n+m,n+m); zeros(n+m,n+m)]; % this is a 18x6 matrix

col2 = [psi; A; zeros(m+n,n)]; % this is a 18x4 matrix

col3 = [-J; zeros(m+n, m+n); FKH]; % this is a 18x6 matrix

aF = [col1, col2, col3]; % this is a 18x16 matrix
aF = aF(1:16,:);
% pad our matrices (the paper is very vague on this point)
aX0 = zeros(16);
aQ = zeros(6);
aX0(7:10,7:10) = X01;
for i = 1:3
%     aX0(i*6-5:i*6-2,i*6-5:i*6-2) = X01;
    aQ(i*2-1:i*2,i*2-1:i*2) = Q;
end
% aF = [0*(FKH + J),        psi,             
% -J;                % this row is 6x6, 8x4, 6x6
%     zeros(n+m,n+m),      A,                zeros(n+m,n+m);  % this row is 6x6, 4x4, 6x6
%     zeros(n+m,n+m)            zeros(n,n), FKH];           % this row is 6x6  4x4  6x6

% The covariance of the robust state estimate is given by the (n+m)×(n+m) upper left
% block of the matrix --> from the paper
% Solve for dPE = aF*PE + PF.' + g*Q*g.'
[T, P] = RiccatiODE(aF.', zeros(16,1), g*aQ*g.',0, aX0, [0 4]); % NO DISTURBANCE
% [T, P] = RiccatiODE(A, zeros(4,1), Gm*Q*Gm.', 0, [X01, zeros(4,2); zeros(2,6)], tv);

% Finding and plotting the covariance, which determines how effective a
% filter matrix is
Pcov = zeros(6, 6, size(T, 1));
for i = 1:length(T)
    Pcov(:, :, i) = P(1:6, 1:6, i);
end
figure;
for i = 1:4
    plot(T,reshape(Pcov(i,i,:), [size(T,1),1]));
    hold on;
end
xlim([2.2, 2.5])
ylim([-10,1E100])
legend('x1','x2','x3','x4');
xlabel('time (sec)');
ylabel('Estimation Error/Covariance');
hold off;
y = zeros(size(T,1),1);
for i = 1:10:1000
    y(i) = norm(Pcov(1:4,:,i));
end
figure;
plot(T,y)
% [T, P] = ode45(RiccatiSimple(Am, Gm*Q*Gm.'), tv, [X02; 0; 0]);

%% Simulation w/ Disturbance
dA = A - Am;
dC = C - Cm;
psi = [dA; dC*(Am + eye(4,4) + dA)]; % 8x4 matrix
% aF is a x matrix
col1 = [0*(FKH + J); zeros(n+m,n+m); zeros(n+m,n+m)]; % this is a 18x6 matrix
col2 = [psi; A; zeros(m+n,n)]; % this is a 18x4 matrix
col3 = [-J; zeros(m+n, m+n); FKH]; % this is a 18x6 matrix
aF = [col1, col2, col3]; % this is a 18x16 matrix
aF = aF(1:16,:);
[T, P] = RiccatiODE(aF.', [zeros(2); ones(2); zeros(12,2)], g*aQ*g.',0, aX0, [0 4]); % YES DISTURBANCE
Pcov = zeros(6, 6, size(T, 1));
for i = 1:length(T)
    Pcov(:, :, i) = P(1:6, 1:6, i);
end
figure;
for i = 1:4
    plot(T,reshape(Pcov(i,i,:), [size(T,1),1]));
    hold on;
end
xlim([2.2, 2.5])
ylim([-10,1E100])
legend('x1','x2','x3','x4');
xlabel('time (sec)');
ylabel('Estimation Error/Covariance');
hold off;
y = zeros(size(T,1),1);
for i = 1:10:1000
    y(i) = norm(Pcov(1:4,:,i));
end
figure;
plot(T,y)
%% Robust Error Dynamics Estimator (failed)
% sys2_1 = ss(Am,B1,Cm,D1);
% figure
% step(sys2_1)
% Udist = randn(4,size(tv,2));
% Unoise = randn(size(tv));
% U = 0*tv;
% augM = ones(1,200);
% [y, t] =  lsim(sys2_1, augM, tv);
% figure(3)
% plot(t, y)
% 
% %% Simulation I Cont'd
% FKH = F-K1.'*H*Zmt;
% 
% % Equation 22
% for i = tv
%     z = ones(6,1)*tv;
%     Zmt = FKH*z - K1*zeros(6,1) + Gm*zeros(6,1); % disturbance 
% 
% end
% 
% % Equation 25
% Error = (F-K1*H) + gam1 - nt -K*dV_dt + G*dB_dt;
% sysK = ss(Am-K1*Cm, [B K1], eye(4), 0*[B K1]);
% figure;
% step(sysK);

%% Local function of the Ricatti ODE
function dPdt = RiccatiSimple(A, Q)
dPdt =@(t,P) A.'*P + P*A + Q; % Set up the Ricatti ODE form
end

%% Function to solve Ricatti ODE
% REQUIRES THE IVP Solver MATLAB TOOLBOX
function [t,P] = RiccatiODE(A,B,Q,R,X0,T)
n = size(B,1); % obtain sizes to use with the IVP toolbox functions
dPdt = @(t,P) A.'*P + P*A + Q;
% (P*B)/R*(B.'*P)
% converts matrix ODE into vector ODE
dydt = mat2vec_ODE(dPdt,n);
% converts initial condition from matrix to vector
X0 = mat2vec_IC(X0);
% solves Riccati ODE
[t,P2] = ode45(dydt,T,X0);
% transforms solution vector into solution matrix
P = vec2mat_sol(P2,n);
t = flipud(t);
end