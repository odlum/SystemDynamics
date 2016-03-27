%System Dynamics 6.1
%McCall, Odlum, Rothberg

%Part A:
%Modifying motor model to include coupling

%Part B:
%Compute and plot step response for: 
%Beta' (rotor angular velocity)
%Motor internal torque

%Time vector
t_max = 1;
t_step = 0.0001;
t = linspace(0,t_max,t_max/t_step);

%Input voltage
ea = 10;       %Input is 10VDC

%Motor constants
Ra = 5;        %Armature resistance, V/A
La = 0.00891;  %Armature inductance, V*s/A
Jr = 0.004;    %Motor Moment of Intertia, N*m*s^2
br = 0.001;    %Motor viscous coefficient, N*m*s
Kb = 0.08;     %Back EMF constant
K = 5;         %Motor torque constant, N*m/A
JL = 0.002;    %Load Moment of Inertia, N*m*s^2
bL = 0.005;    %Viscous coefficient of load, N*m*s
ktc= 200;      %Shaft coupler stiffness, N*m

%Rewriting state space equations with coupler
%Electrical system of DC Motor:
    %Ra*ia + La*ia' + kb*theta' = ea
%Mechanical system of DC Motor:
    %Jr*theta" + br*theta'  - K*ia + ktc(theta - beta) = 0
%Coupling between electrical and mechanical motor parts
    %eb = Kb*theta'
    %Tm = K*ia
%EOM of mechanical load:
    %JL*beta" + bL*beta' + ktc(beta - theta) = 0 
    

%State Space Equations:
%x1 = ia
%x2 = theta
%x3 = theta'
%x4 = beta
%x5 = beta'

%input
A = [-Ra/La     0   -Kb/La      0       0;...
        0       0       1       0       0;...
    K/Jr    -ktc/Jr  -br/Jr     ktc/Jr  0;...
        0       0       0       0       1;...
        0   ktc/JL      0   -ktc/JL   -bL/JL];    
B = [1/La; 0; 0; 0; 0];

%output: solving for beta' and Tm
C = [0 0 0 0 1; K 0 0 0 0];
D = [0; 0];

%perform state space operation
sys = ss(A,B,C,D);

%define input voltage step function
y = ea*step(sys,t);

%plot responses
figure(1)
plot(t,y(:,1))
    title('Response of $\dot{\beta}$ vs. Time','interpreter','latex')
    xlabel('Time, s')
    ylabel('\omega, rad/s')
    grid
    
figure(2)
plot(t,y(:,2))
    title('Response of Internal Torque vs. Time')
    xlabel('Time, s')
    ylabel('T_m, Nm')
    grid




