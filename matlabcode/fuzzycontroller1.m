function [sys,x0,str,ts] = spacemodel(t,x,u,flag)
switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1,
    sys=mdlDerivatives(t,x,u);
case 3,
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 2*9;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 12;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = zeros(1,2*9);
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)
q1_desired = u(1); % desired angular position of joint 1
dq1_desired = u(2); % desired angular velocity of joint 1
ddq1_desired = u(3); % desired angular acceleration of joint 1
q2_desired = u(4); % desired angular position of joint 2
dq2_desired = u(5); % desired angular velocity of joint 2
ddq2_desired = u(6); % desired angular acceleration of joint 2

q1 = u(7); % actual angular position of joint 1
dq1 = u(8); % actual angular velocity of joint 1
q2 = u(9); % actual angular position of joint 2
dq2 = u(10); % actual angular velocity of joint 2

dq1_delay = u(11); % actual angular velocity of joint 1 at last sampling moment
dq2_delay = u(12); % actual angular velocity of joint 2 at last sampling moment

Ts = 0.001; % sampling period
ddq1 = (dq1 - dq1_delay)/Ts; % actual angular acceleration of joint 1
ddq2 = (dq2 - dq2_delay)/Ts; % actual angular acceleration of joint 2

error1 = q1_desired - q1; % tracking error of joint 1
error2 = q2_desired - q2; % tracking error of joint 2
error1_derivative = dq1_desired - dq1; % derivative of tracking error of joint 1
error2_derivative = dq2_desired - dq2; % derivative of tracking error of joint 2
error1_second_derivative = ddq1_desired - ddq1; % second order derivative of tracking error of joint 1
error2_second_derivative = ddq2_desired - ddq2; % second order derivative of tracking error of joint 2

z1 = [error1; error1_derivative]; % first fuzzy controller input
z2 = [error2; error2_derivative]; % second fuzzy controller input

% membership function of first fuzzy controller
membership1 = struct();

% membership function of error1
membership1.error1.negative = exp(-1/2 * ((z1(1) + 1.25)/0.6)^2);
membership1.error1.medium = exp(-1/2 * ((z1(1))/0.6)^2);
membership1.error1.positive = exp(-1/2 * ((z1(1) - 1.25)/0.6)^2);

% membership function of error1_derivative
membership1.error1_derivative.negative = exp(-1/2 * ((z1(2) + 1.25)/0.6)^2);
membership1.error1_derivative.medium = exp(-1/2 * ((z1(2))/0.6)^2);
membership1.error1_derivative.positive = exp(-1/2 * ((z1(2) - 1.25)/0.6)^2);

% membership function of second fuzzy controller
membership2 = struct();

% membership function of error2
membership2.error2.negative = exp(-1/2 * ((z2(1) + 1.25)/0.6)^2);
membership2.error2.medium = exp(-1/2 * ((z2(1))/0.6)^2);
membership2.error2.positive = exp(-1/2 * ((z2(1) - 1.25)/0.6)^2);

% membership function of error2_derivative
membership2.error2_derivative.negative = exp(-1/2 * ((z2(1) + 1.25)/0.6)^2);
membership2.error2_derivative.medium = exp(-1/2 * ((z2(1))/0.6)^2);
membership2.error2_derivative.positive = exp(-1/2 * ((z1(1) - 1.25)/0.6)^2);

% fuzzy basis function xi_1 of first fuzzy controller
denominator_xi1 = 0;
numerator_xi1 = zeros(1,9);
for i=1:1:3
    for j=1:1:3
        fields1 = fieldnames(membership1.error1);
        fieldName1 = fields1{i};
        mu1 = membership1.error1.(fieldName1);
        
        fields2 = fieldnames(membership1.error1_derivative);
        fieldName2 = fields2{j};
        mu2 = membership1.error1_derivative.(fieldName2);
        
        numerator_xi1(3*(i-1)+j) = mu1 * mu2; % numerator of fuzzy basis function xi_1
        denominator_xi1 = denominator_xi1 + mu1 * mu2; % denominator of fuzzy basis function xi_1
    end
end

xi_1 = numerator_xi1/denominator_xi1;

% fuzzy basis function xi_1 of second fuzzy controller
denominator_xi2 = 0;
numerator_xi2 = zeros(1,9);
for i=1:1:3
    for j=1:1:3
        fields1 = fieldnames(membership2.error2);
        fieldName1 = fields1{i};
        mu1 = membership2.error2.(fieldName1);
        
        fields2 = fieldnames(membership2.error2_derivative);
        fieldName2 = fields2{j};
        mu2 = membership2.error2_derivative.(fieldName2);
        
        numerator_xi2(3*(i-1)+j) = mu1 * mu2; % numerator of fuzzy basis function xi_2
        denominator_xi2 = denominator_xi2 + mu1 * mu2; % denominator of fuzzy basis function xi_2
    end
end

xi_2 = numerator_xi2/denominator_xi2;

lambda1 = 1.0;
lambda2 = 1.0;

error1_state = error1 + lambda1 * error1_derivative; % filtered tracking error of joint 1
error2_state = error2 + lambda2 * error2_derivative; % filtered tracking error of joint 2
error1_state_derivative = error1_derivative + lambda1 * error1_second_derivative; % derivative of filtered tracking error of joint 1
error2_state_derivative = error2_derivative + lambda2 * error2_second_derivative; % derivative of filtered tracking error of joint 2

eta = 10; % learning rate
epsilon0 = 0.01;
K = [20, 0;0, 10]; % control gain
K0 = [20, 0;0, 10];
delta = 0.001;

theta1 = zeros(9, 1);
theta2 = zeros(9, 1);
for i = 1:9
    theta1(i) = x(i); % first fuzzy controller parameter
    theta2(i) = x(i + 9);% second fuzzy controller parameter
end

theta1 = theta1';
theta2 = theta2';

% derivative of fuzzy controller parameter
theta1_derivative = eta * xi_1 * (error1_state_derivative + K(1,1) * error1_state + K0(1,1) * tanh(error1_state./epsilon0)) - eta * delta * theta1;
theta2_derivative = eta * xi_2 * (error2_state_derivative + K(2,2) * error2_state + K0(2,2) * tanh(error2_state./epsilon0)) - eta * delta * theta2;

for i = 1:9
    sys(i) = theta1_derivative(i);
end

for i = 1:9
    sys(i+9) = theta2_derivative(i);
end

function sys=mdlOutputs(t,x,u)
q1_desired = u(1); % desired angular position of joint 1
dq1_desired = u(2); % desired angular velocity of joint 1
q2_desired = u(4); % desired angular position of joint 2
dq2_desired = u(5); % desired angular velocity of joint 2

q1 = u(7); % actual angular position of joint 1
dq1 = u(8); % actual angular velocity of joint 1
q2 = u(9); % actual angular position of joint 2
dq2 = u(10); % actual angular velocity of joint 2

error1 = q1_desired - q1; % tracking error
error2 = q2_desired - q2;
error1_derivative = dq1_desired - dq1; % derivative of tracking error
error2_derivative = dq2_desired - dq2;

% fuzzy controller input
z1 = [error1; error1_derivative];
z2 = [error2; error2_derivative];

% membership function of first fuzzy controller
membership1 = struct();

% membership function of error1
membership1.error1.negative = exp(-1/2 * ((z1(1) + 1.25)/0.6)^2);
membership1.error1.medium = exp(-1/2 * ((z1(1))/0.6)^2);
membership1.error1.positive = exp(-1/2 * ((z1(1) - 1.25)/0.6)^2);

% membership function of error1_derivative
membership1.error1_derivative.negative = exp(-1/2 * ((z1(2) + 1.25)/0.6)^2);
membership1.error1_derivative.medium = exp(-1/2 * ((z1(2))/0.6)^2);
membership1.error1_derivative.positive = exp(-1/2 * ((z1(2) - 1.25)/0.6)^2);

% membership function of second fuzzy controller
membership2 = struct();

% membership function of error2
membership2.error2.negative = exp(-1/2 * ((z2(1) + 1.25)/0.6)^2);
membership2.error2.medium = exp(-1/2 * ((z2(1))/0.6)^2);
membership2.error2.positive = exp(-1/2 * ((z2(1) - 1.25)/0.6)^2);

% membership function of error2_derivative
membership2.error2_derivative.negative = exp(-1/2 * ((z2(1) + 1.25)/0.6)^2);
membership2.error2_derivative.medium = exp(-1/2 * ((z2(1))/0.6)^2);
membership2.error2_derivative.positive = exp(-1/2 * ((z1(1) - 1.25)/0.6)^2);

% fuzzy basis function xi_1 of first fuzzy controller
denominator_xi1 = 0;
numerator_xi1 = zeros(1,9);
for i=1:1:3
    for j=1:1:3
        fields1 = fieldnames(membership1.error1);
        fieldName1 = fields1{i};
        mu1 = membership1.error1.(fieldName1);
        
        fields2 = fieldnames(membership1.error1_derivative);
        fieldName2 = fields2{j};
        mu2 = membership1.error1_derivative.(fieldName2);
        
        numerator_xi1(3*(i-1)+j) = mu1 * mu2; % numerator of fuzzy basis function xi_1
        denominator_xi1 = denominator_xi1 + mu1 * mu2; % denominator of fuzzy basis function xi_1
    end
end

xi_1 = numerator_xi1/denominator_xi1;

% fuzzy basis function xi_1 of second fuzzy controller
denominator_xi2 = 0;
numerator_xi2 = zeros(1,9);
for i=1:1:3
    for j=1:1:3
        fields1 = fieldnames(membership2.error2);
        fieldName1 = fields1{i};
        mu1 = membership2.error2.(fieldName1);
        
        fields2 = fieldnames(membership2.error2_derivative);
        fieldName2 = fields2{j};
        mu2 = membership2.error2_derivative.(fieldName2);
        
        numerator_xi2(3*(i-1)+j) = mu1 * mu2; % numerator of fuzzy basis function xi_2
        denominator_xi2 = denominator_xi2 + mu1 * mu2; % denominator of fuzzy basis function xi_2
    end
end

xi_2 = numerator_xi2/denominator_xi2;

theta1 = zeros(9, 1);
theta2 = zeros(9, 1);
for i = 1:9
    theta1(i) = x(i); % first fuzzy controller parameter
    theta2(i) = x(i + 9);% second fuzzy controller parameter
end

% fuzzy controller output, Eq.16
% control amount is equal to fuzzy basis function xi multiplied by fuzzy controller parameter theta
control1 = xi_1 * theta1;
control2 = xi_2 * theta2;

sys(1) = control1;
sys(2) = control2;