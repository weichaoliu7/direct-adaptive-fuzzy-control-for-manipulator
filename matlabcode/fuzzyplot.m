close all;

figure(1);
subplot(211);
plot(t,q1(:,1),'r',t,q1(:,2),'b');
xlabel('time(s)');ylabel('position tracking of joint 1');
h=legend('desired position','reference model position');
set(h,'Box','off');
subplot(212);
plot(t,q2(:,1),'r',t,q2(:,2),'b');
xlabel('time(s)');ylabel('position tracking of joint 2');
h=legend('desired position','actual position');
set(h,'Box','off');
sgtitle('position tracking in joint space');

figure(2);
subplot(211);
plot(t,dq1(:,1),'r',t,dq1(:,2),'b');
xlabel('time(s)');ylabel('velocity tracking of joint 1');
h=legend('desired velocity','actual velocity');
set(h,'Box','off');
subplot(212);
plot(t,dq2(:,1),'r',t,dq2(:,2),'b');
xlabel('time(s)');ylabel('velocity tracking of joint 2');
h=legend('desired velocity','actual velocity');
set(h,'Box','off');
sgtitle('velocity tracking in joint space');

figure(3);
subplot(211);
plot(t,q1(:,2)-q1(:,1),'b');
xlabel('time(s)');ylabel('position error of joint 1');
subplot(212);
plot(t,q2(:,2)-q2(:,1),'b');
xlabel('time(s)');ylabel('position error of joint 2');
sgtitle('position error');

figure(4);
subplot(211);
plot(t,dq1(:,2)-dq1(:,1),'b');
xlabel('time(s)');ylabel('velocity error of joint 1');
subplot(212);
plot(t,dq2(:,2)-dq2(:,1),'b');
xlabel('time(s)');ylabel('velocity error of joint 2');
sgtitle('velocity error');

figure(5);
subplot(211);
plot(t,control(:,1),'b');
xlabel('time(s)');ylabel('torque1');
subplot(212);
plot(t,control(:,2),'b');
xlabel('time(s)');ylabel('torque2');
sgtitle('control input torque');
