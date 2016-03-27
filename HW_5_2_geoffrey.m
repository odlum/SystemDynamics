%System Dynamics HW 5.2
%State Space Calculations
%McCall, Odlum, Rothberg

%define variables
k = 10000; % N/m
b = 5; %Ns/m
m = 10; %kg

t_max = 15; %15 sec time interval
t_step = 0.001;  %time step for graphing
t = linspace(0,t_max,t_max / t_step);

%define state space matrices
%y1 = x1; y2 = y1" = x2'

% %state matrix (n x n)    
A = [0 1 0 0 ;...
     -2*k/m -2*b/m k/m b/m ;...
     0 0 0 1 ;...
     k/m b/m -2*k/m -b/m ];

% %input matrix (n x r)
B = [0 0;...                           
    (1/m) 0;...
    0 0;...
    0 (1/m)]; 

% %output matrix (m x n)
C = [1 0 0 0;...
    0 1 0 0;...
    (-2*k/m) (-2*b/m) (k/m) (b/m)];

% %state direct transmission matrix (m x r)    
D = [0 0; 0 0; (1/m) 0];                                 

%--------------------------------------------------------------------------

%Part A: Plot output for:
%   f1(t) = f2(t) = 0; 
%   y1(0) = 0.1m; 
%   y1'(0) = y2(0) = y2'(0) = 0;

%X0 parameter is initial conditions for 'initial' function
x0 = [0.1*m 0 0 0];

[y,x,t] = initial(A, B, C, D, x0, t);

%Output 1
y1 = y(:,1); 

%Output 2
y2 = y(:,2); 

%Output 3
y3 = y(:,3); 

%--------------------------------------------------------------------------

%Part B: 
%f1(t) = 100 x {(1(t) - 1(t-2)}
f1 = 0*t;
t_pulse = 2 / t_step;
for idx = (1:t_pulse)    
    f1(idx) = 100;
end

sys = ss(A, B, C, D);

%f2(t) = 20 x d(t-3)
f2 = 0*t;
t_delta = 3 / t_step;
f_delta = 20*impulse(sys,t);
% for idx = (3/t_step):(t_max/t_step)    
%     f2(idx) = f_delta(idx - 3/t_step + 1);
% end
f2(t_delta) = 20/t_step;


%y1(0) = y2(0) = y3(0) = y4(0) = 0
x0_b = [0 0 0 0];

u = [f1; f2];

[z,t] = lsim(sys, u, t, x0_b);

%Output 1
z1 = z(:,1); 

%Output 2
z2 = z(:,2); 

%Output 3
z3 = z(:,3); 

%--------------------------------------------------------------------------

% Plot responses
figure (1)

subplot(3,1,1)
plot(t,y1)
ylabel('y1(t)')
grid

subplot(3,1,2)
plot(t,y2)
ylabel('y1"(t)')
grid

subplot(3,1,3)
plot(t,y3)
ylabel('y1""(t)')
xlabel('time (s)')
grid

figure (2)

subplot(3,1,1)
plot(t,z1)
ylabel('y1(t)')
grid

subplot(3,1,2)
plot(t,z2)
ylabel('y1"(t)')
grid

subplot(3,1,3)
plot(t,z3)
ylabel('y1""(t)')
xlabel('time (s)')
grid
