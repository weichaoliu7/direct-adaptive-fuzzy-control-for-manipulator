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
control = [u(1);u(2)];

m1 = 1;
me = 2;
l1 = 1;
lc1 = 0.5;
lce = 0.6;
I1 = 0.12;
Ie = 0.25;
theta_e = pi/6;

a1 = I1 + m1 * lc1^2 + Ie + me * lce^2 + me * l1^2;
a2 = Ie + me * lce^2;
a3 = me * l1 *lce * cos(theta_e);
a4 = me * l1 * lce * sin(theta_e);

% inertia matrix for manipulator dynamics equation
M(1,1) = a1 + 2 * a3 * cos(x(3)) + 2 * a4 * sin(x(3));
M(1,2) = a2 + a3 * cos(x(3)) + a4 * sin(x(3));
M(2,1) = M(1,2);
M(2,2) = a2;

% Coriolis and centrifugal matrix for manipulator dynamics equation
h = a3 * sin(x(3)) - a4 * cos(x(3));
C(1,1) = -h * x(4);
C(1,2) = -h * (x(2) + x(4));
C(2,1) = h * x(2);
C(2,2) = 0;

dq = [x(2);x(4)]; % actual angular velocity of manipulator
f = - inv(M) * C * dq; % unknown smooth nonlinear function f(x)
G = inv(M); % unknown smooth nonlinear function G(x)
S = f + G * control; % angular acceleration, manipulator dynamic system compact form, Eq. 2

sys(1)=x(2); % dq1
sys(2)=S(1); % ddq1
sys(3)=x(4); % dq2
sys(4)=S(2); % ddq2
function sys=mdlOutputs(t,x,u)
sys(1)=x(1); % q1
sys(2)=x(2); % dq1
sys(3)=x(3); % q2
sys(4)=x(4); % dq2