function [sys,x0,str,ts] = spacemodel(t,x,u,flag)
switch flag
case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1
    sys=mdlDerivatives(t,x,u);
case 3
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 1;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 12;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = 0;
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)
sys=0;

function sys=mdlOutputs(t,x,u)
q1_desired = u(1); % desired angular position of joint 1
dq1_desired = u(2); % desired angular velocity of joint 1
q2_desired = u(4); % desired angular position of joint 2
dq2_desired = u(5); % desired angular velocity of joint 2

q1 = u(7); % actual angular position of joint 1
dq1 = u(8); % actual angular velocity of joint 1
q2 = u(9); % actual angular position of joint 2
dq2 = u(10); % actual angular velocity of joint 2

error1 = q1_desired - q1; % position tracking error of joint 1
error2 = q2_desired - q2; % position tracking error of joint 2
error1_derivative = dq1_desired - dq1; % velocity tracking error of joint 1
error2_derivative = dq2_desired - dq2; % velocity tracking error of joint 2
error = [error1;error2];
error_derivative = [error1_derivative;error2_derivative];
Kd = [100,0;0,100];
Kp = [2000,0;0,2000];

% control input torque
torque = Kp*error + Kd*error_derivative; % PD controller

sys(1)=torque(1);
sys(2)=torque(2);