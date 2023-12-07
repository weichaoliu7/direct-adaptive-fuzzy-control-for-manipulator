function [sys,x0,str,ts] = s_function(t,x,u,flag)

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
sizes.NumContStates  = 4;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 4;
sizes.NumInputs      = 2;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;
sys=simsizes(sizes);
x0=[0.0 2.0 0.0 2.5]; % manipulator initial condition
str=[];
ts=[];

function sys=mdlDerivatives(t,x,u)
l = 1; % 1ength of link
m1 = 1; % mass of link 1
m2 = 2; % mass of link 2
g = 9.41; % gravitational acceleration

kv1=0.3; % coefficient of viscous friction
kc1=0.2; % coefficient of coulomb friction
kv2=0.5;
kc2=0.5;

% inertia matrix for manipulator dynamics equation
D(1,1) = 1/3 * m1 * l^2 + 4/3 * m2 * l^2 + m2 * l^2 * cos(x(3));
D(1,2) = 1/3 * m2 * l^2 + 1/2 * m2 * l^2 * cos(x(3));
D(2,1) = D(1,2);
D(2,2) = 1/3 * m2 * l^2;

% Coriolis and centrifugal matrix for manipulator dynamics equation 
C(1,1) = -m2 * l^2 * sin(x(3)) * x(4);
C(1,2) = -1/2 * m2 * l^2 * sin(x(3)) * x(4);
C(2,1) = 1/2 * m2 * l^2 * sin(x(3)) * x(2);
C(2,2) = 0;

% gravitational torque matrix for manipulator dynamics equation 
G(1) = 1/2 * m1 * g * l * cos(x(1)) + 1/2 * m2 * g * l * cos(x(1)+x(3)) + m2 * g * l * cos(x(1));
G(2) = 1/2 * m2 * g * l * cos(x(1)+x(3));

Fv(1) = kv1 * x(2); % viscous friction torque
Fv(2) = kv2 * x(4);

Fc(1) = kc1 * sign(x(2)); % coulomb friction torque
Fc(2) = kc2 * sign(x(4));

torque=u(1:2); % control input torque
dq=[x(2);x(4)];
S=inv(D)*(torque-C*dq-G'-Fv'-Fc'); % angular acceleration, from manipulator dynamics equation

sys(1)=x(2); % dq1
sys(2)=S(1); % ddq1
sys(3)=x(4); % dq2
sys(4)=S(2); % ddq2
function sys=mdlOutputs(t,x,u)
sys(1)=x(1); % q1
sys(2)=x(2); % dq1
sys(3)=x(3); % q2
sys(4)=x(4); % dq2