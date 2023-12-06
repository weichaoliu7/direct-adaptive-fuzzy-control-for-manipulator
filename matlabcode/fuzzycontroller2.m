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
sizes.NumContStates  = 2*3^4;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 12;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0  = zeros(1,2*3^4);
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

z = [error1;error1_derivative;error2;error2_derivative]; % fuzzy controller input

% membership function of fuzzy controller
membership = struct();

% membership function of error1
membership.error1.negative = exp(-1/2 * ((z(1) + 1.25)/0.6)^2);
membership.error1.medium = exp(-1/2 * ((z(1))/0.6)^2);
membership.error1.positive = exp(-1/2 * ((z(1) - 1.25)/0.6)^2);

% membership function of error1_derivative
membership.error1_derivative.negative = exp(-1/2 * ((z(2) + 1.25)/0.6)^2);
membership.error1_derivative.medium = exp(-1/2 * ((z(2))/0.6)^2);
membership.error1_derivative.positive = exp(-1/2 * ((z(2) - 1.25)/0.6)^2);

% membership function of error2
membership.error2.negative = exp(-1/2 * ((z(3) + 1.25)/0.6)^2);
membership.error2.medium = exp(-1/2 * ((z(3))/0.6)^2);
membership.error2.positive = exp(-1/2 * ((z(3) - 1.25)/0.6)^2);

% membership function of error2_derivative
membership.error2_derivative.negative = exp(-1/2 * ((z(4) + 1.25)/0.6)^2);
membership.error2_derivative.medium = exp(-1/2 * ((z(4))/0.6)^2);
membership.error2_derivative.positive = exp(-1/2 * ((z(4) - 1.25)/0.6)^2);

% fuzzy basis function xi of fuzzy controller
denominator_xi = 0;
numerator_xi = zeros(1,9);
for i=1:1:3
    for j=1:1:3
        for k=1:1:3
            for l=1:1:3
                fields1 = fieldnames(membership.error1);
                fieldName1 = fields1{i};
                mu1 = membership.error1.(fieldName1);

                fields2 = fieldnames(membership.error1_derivative);
                fieldName2 = fields2{j};
                mu2 = membership.error1_derivative.(fieldName2);

                fields3 = fieldnames(membership.error2);
                fieldName3 = fields3{k};
                mu3 = membership.error2.(fieldName3);

                fields4 = fieldnames(membership.error2_derivative);
                fieldName4 = fields4{l};
                mu4 = membership.error2_derivative.(fieldName4);

                numerator_xi(3 * (( 3 * (( 3 * (i-1) + j ) -1 ) + k) - 1) + l) = mu1 * mu2 * mu3 * mu4; % numerator of fuzzy basis function xi
                denominator_xi = denominator_xi + mu1 * mu2 * mu3 * mu4; % denominator of fuzzy basis function xi
            end
        end
    end
end

xi = numerator_xi/denominator_xi;

lambda1 = 1.0;
lambda2 = 1.0;

error1_state = error1 + lambda1 * error1_derivative; % filtered tracking error of joint 1
error2_state = error2 + lambda2 * error2_derivative; % filtered tracking error of joint 2
error_state = [error1_state; error2_state]; % filtered tracking error

error1_state_derivative = error1_derivative + lambda1 * error1_second_derivative; % derivative of filtered tracking error of joint 1
error2_state_derivative = error2_derivative + lambda2 * error2_second_derivative; % derivative of filtered tracking error of joint 2
error_state_derivative = [error1_state_derivative;error2_state_derivative]; % derivative of filtered tracking error

eta = 10; % learning rate
epsilon0 = 0.01; % small positive constant
K = [100, 0;0, 50]; % control gain
K0 = [100, 0;0, 50];
delta = 0.001; % delta-modification term parameter

theta = zeros(81,2);
for i = 1:81
    theta(i,1) = x(i); % fuzzy controller parameter
    theta(i,2) = x(i + 81);
end

% derivative of fuzzy controller parameter
theta_derivative = eta * xi' * (error_state_derivative + K * error_state + K0 * tanh(error_state./epsilon0))' - eta * delta * theta;

for i = 1:81
    sys(i) = theta_derivative(i,1);
    sys(i+81) = theta_derivative(i,2);
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

z = [error1;error1_derivative;error2;error2_derivative]; % fuzzy controller input

% membership function of fuzzy controller
membership = struct();

% membership function of error1
membership.error1.negative = exp(-1/2 * ((z(1) + 1.25)/0.6)^2);
membership.error1.medium = exp(-1/2 * ((z(1))/0.6)^2);
membership.error1.positive = exp(-1/2 * ((z(1) - 1.25)/0.6)^2);

% membership function of error1_derivative
membership.error1_derivative.negative = exp(-1/2 * ((z(2) + 1.25)/0.6)^2);
membership.error1_derivative.medium = exp(-1/2 * ((z(2))/0.6)^2);
membership.error1_derivative.positive = exp(-1/2 * ((z(2) - 1.25)/0.6)^2);

% membership function of error2
membership.error2.negative = exp(-1/2 * ((z(3) + 1.25)/0.6)^2);
membership.error2.medium = exp(-1/2 * ((z(3))/0.6)^2);
membership.error2.positive = exp(-1/2 * ((z(3) - 1.25)/0.6)^2);

% membership function of error2_derivative
membership.error2_derivative.negative = exp(-1/2 * ((z(4) + 1.25)/0.6)^2);
membership.error2_derivative.medium = exp(-1/2 * ((z(4))/0.6)^2);
membership.error2_derivative.positive = exp(-1/2 * ((z(4) - 1.25)/0.6)^2);

% fuzzy basis function xi of fuzzy controller
denominator_xi = 0;
numerator_xi = zeros(1,9);
for i=1:1:3
    for j=1:1:3
        for k=1:1:3
            for l=1:1:3
                fields1 = fieldnames(membership.error1);
                fieldName1 = fields1{i};
                mu1 = membership.error1.(fieldName1);

                fields2 = fieldnames(membership.error1_derivative);
                fieldName2 = fields2{j};
                mu2 = membership.error1_derivative.(fieldName2);

                fields3 = fieldnames(membership.error2);
                fieldName3 = fields3{k};
                mu3 = membership.error2.(fieldName3);

                fields4 = fieldnames(membership.error2_derivative);
                fieldName4 = fields4{l};
                mu4 = membership.error2_derivative.(fieldName4);

                numerator_xi(3 * (( 3 * (( 3 * (i-1) + j ) -1 ) + k) - 1) + l) = mu1 * mu2 * mu3 * mu4; % numerator of fuzzy basis function xi
                denominator_xi = denominator_xi + mu1 * mu2 * mu3 * mu4; % denominator of fuzzy basis function xi
            end
        end
    end
end

xi = numerator_xi/denominator_xi;

theta = zeros(81,2);
for i = 1:81
    theta(i,1) = x(i); % fuzzy controller parameter
    theta(i,2) = x(i + 81);
end

% fuzzy controller output, Eq.16
% control amount is equal to fuzzy basis function xi multiplied by fuzzy controller parameter theta
control = xi * theta;

sys(1) = control(1);
sys(2) = control(2);