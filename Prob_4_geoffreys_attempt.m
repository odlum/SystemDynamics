%HW 2.4
%Transient Response Analysis with MATLAB
%Odlum, Geoffrey; Rothberg, Benjamin; McCall, Troy
%
close all
clear all
clc

%--------------------------------------------------------------------------

%Define System Variables
m = 10;                         %Mass
c = 5;                          %Damping Constant
k = 50;                         %Spring Constant

freq_nat = sqrt(k/m);           %Natural Frequency
T_nat = 2*pi/freq_nat;          %Natural Period

imp = 4 * 0.75;                 %Magnitude of Impulse
t_step = 0.01;                  %step for time vectors
t_max = 12;                     %we want response from 0 to 12s

%--------------------------------------------------------------------------

%PART 1 : response x(t)
%Define System in Terms of Transfer Function
num = [1];
den = [m,c,k];
sys = tf(num, den);

%Response of System
t_pulse = 0.75;                 %Duration of Pulse
t = [0:t_step:t_max];           %Create time vector
f = 0*t;                        %Create Zero Vector from t

%Populate Zero Vector with impulse divided by pulse duration
for idx = 1:t_pulse/t_step
    f(idx) = imp/t_pulse;       
end

x_pulse = lsim(sys,f,t);        %Compute Response of System

%Compute Response of System Using Delta Funciton
d_imp = 2;                          %delta impulse magnitude
dt = 1.2;                           %delta impulse time
x_delta = d_imp*impulse(sys,t);     %delta function
x_delta_1 = 0*x_delta;              %vector for shifted delta function

%loop populates x_delta_1 with values from x_delta, shifted 1.2s
for idx = (dt/t_step):(t_max/t_step)    
    x_delta_1(idx) = x_delta(idx - dt/t_step + 1);
end

%total displacement response is sum of individual responses
x_tot = x_delta_1 + x_pulse;

%--------------------------------------------------------------------------

%PART 2 : response x'(t)
%compute the velocity response of the system
%Define System in Terms of Transfer Function
num_1 = [1,0];
den_1 = [m,c,k];
sys_1 = tf(num_1, den_1);

v_pulse = lsim(sys_1,f,t);        %Compute Response of System

%Compute Response of System Using Delta Funciton
v_delta = d_imp*impulse(sys_1,t);   %delta function
v_delta_1 = 0*v_delta;              %vector for shifted delta function

%for loop to populate v_delta_1 vector with values shifted 1.2s
for idx = (dt/t_step):(t_max/t_step)
    v_delta_1(idx) = v_delta(idx - dt/t_step + 1);
end

%total velocity response is sum of individual responses
v_tot = v_delta_1 + v_pulse;

%--------------------------------------------------------------------------

%PART B : suppressing responses x(t) and x'(t) due to initial pulse
%first step is to find one of the zeros of response x(t)

%loop populates x_delta_1 with values from x_delta, shifted 1.2s
v_pulse_1 = v_pulse;
for idx = 1:(t_pulse/t_step)    
    v_pulse_1(idx) = 0;
end

%find the first instance after t = 0.75s where x(t) = 0

%while loop searches through x_pulse response until x drops below 0
idx = 2;
while x_pulse(idx) > 0
    idx = idx + 1;
end

%suppressing time is set
t_min = idx - 1;

%suppressing velocity is determined
v_max = v_pulse(t_min);

%determine magnitude of momentum at this point
%MV = FT = I; 
%imp_sup is the magnitude of the impulse needed to suppress motion
imp_sup = -v_max * m;

%--------------------------------------------------------------------------

%Compute Response of System Using Delta Function
%x(t) response
x_sup = imp_sup*impulse(sys,t);       %delta function
x_sup_1 = 0*x_sup;              %vector for shifted delta function

%loop populates x_sup_1 with values from x_sup, shifted t_min s
for idx = (t_min):(t_max/t_step)    
    x_sup_1(idx) = x_sup(idx - t_min + 1);
end

%total displacement response is sum of individual responses
x_tot_sup = x_sup_1 + x_pulse;

%x'(t) response
%Compute Response of System Using Delta Funciton
v_sup = imp_sup*impulse(sys_1,t);   %delta function
v_sup_1 = 0*v_delta;              %vector for shifted delta function

%for loop to populate v_delta_1 vector with values shifted 1.2s
for idx = (t_min):(t_max/t_step)
    v_sup_1(idx) = v_sup(idx - t_min + 1);
end

%total velocity response is sum of individual responses
v_tot_sup = v_sup_1 + v_pulse;
%--------------------------------------------------------------------------

%CONCLUSION : graphing A and B
%Graph Responses
figure(1)   %part A
plot(t,x_tot,t,v_tot)
    grid;
    title('Response of x(t), x''(t) to f1 + f2');
    xlabel('Time (seconds)');
    ylabel('x(t)');
    legend('x(t)', 'x''(t)');

figure(2)   %part B
plot(t,x_tot_sup,t,v_tot_sup)
    grid;
    title('Response of x(t), x''(t) to f1 + suppressing delta function');
    xlabel('Time (seconds)');
    ylabel('x(t)');
    legend('x(t) suppressed', 'x''(t) suppressed');
    
    