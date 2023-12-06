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
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [0.001 0];
function sys=mdlOutputs(t,x,u)
q1_desired=pi*(sin(0.5*t)+0.1*sin(2*t)); % desired angular position of joint 1
dq1_desired=pi*(0.5*cos(0.5*t)+0.1*2*cos(2*t)); % desired angular velocity of joint 1
ddq1_desired=-pi*((0.5^2)*sin(0.5*t)+0.1*(2^2)*sin(2*t)); % desired angular acceleration of joint 1
q2_desired=pi*(0.5*sin(t)+0.1*sin(3*t)); % desired angular position of joint 2
dq2_desired=pi*(0.5*cos(t)+0.1*3*cos(3*t)); % desired angular velocity of joint 2
ddq2_desired=-pi*(0.5*sin(t)+0.1*(3^2)*sin(3*t)); % desired angular acceleration of joint 2

sys(1)=q1_desired;
sys(2)=dq1_desired;
sys(3)=ddq1_desired;
sys(4)=q2_desired;
sys(5)=dq2_desired;
sys(6)=ddq2_desired;