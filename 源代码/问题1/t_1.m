clear,clc
a=readmatrix('../../赛题/附件3.xlsx');
a=a(:,2:end);

i = 1;
w = a(i, 1);
m_a1 = a(i, 2);
c1 = a(i, 4);
f = a(i, 6);

m_f = 4866;
r = 1;
m_o = 2433;
rho = 1025;
k = 80000;
g = 9.8;
l0 = 0.5;

xf0 = -2;
xo0 = -1.8;

c = 10000;

T = 2 * pi / w;

dp_c=@(t,p,c)[p(2);
          (-(c1 + c) * p(2) - k * p(1) + c * p(4) + k * p(3) - k * l0 + ...
          rho * g * V(p(1)) - m_f * g +...
          f * cos(w * t)) / (m_a1 + m_f);
          p(4);
          (c * p(4) + k * p(3) - c * p(2) - k * p(1) + m_o * g - k * l0) / (-m_o);];
dp1 = @(t,p)dp_c(t, p, c);
dp2 = @(t,p)dp_c(t, p, c*(abs(p(2) - p(4))^0.5));

sol1=ode45(dp1,[0,40*T],[-2 0 -1.8 0]);
sol2=ode45(dp2,[0,40*T],[-2 0 -1.8 0]);


%保存结果
t = 0:0.2:40*T;
t_dis = [10, 20, 40, 60, 100];
% result1-1
p1 = deval(sol1,t);
xf1 = p1(1,:);
vf1 = p1(2,:);
xo1 = p1(3,:);
vo1 = p1(4,:);

figure;
plot(t, xf1);
xlabel('$t/s$', 'Interpreter', 'latex', 'fontsize',16)
ylabel('$x_f/m$', 'Interpreter', 'latex', 'fontsize',16)

figure;
plot(t, xo1, 'r');
xlabel('$t/s$', 'Interpreter', 'latex', 'fontsize',16)
ylabel('$x_o/m$', 'Interpreter', 'latex', 'fontsize',16)

result1 = [t; xf1; vf1; xo1; vo1]';
filename = '../../结果/result1-1.xlsx';
writematrix(result1,filename,'Sheet',1,'Range','A3:E900')

% 论文展现
p1_dis = deval(sol1,t_dis);
xf1_dis = p1_dis(1,:);
vf1_dis = p1_dis(2,:);
xo1_dis = p1_dis(3,:);
vo1_dis = p1_dis(4,:);
result1_dis = [t_dis; xf1_dis; vf1_dis; xo1_dis; vo1_dis]';
filename = '../../结果/result1-1_dis.xlsx';
writematrix(result1_dis,filename,'Sheet',1,'Range','A3:E7')

% result1-2
p2 = deval(sol2,t);
xf2 = p2(1,:);
vf2 = p2(2,:);
xo2 = p2(3,:);
vo2 = p2(4,:);

figure;
plot(t, xf2);
xlabel('$t/s$', 'Interpreter', 'latex', 'fontsize',16)
ylabel('$x_f/m$', 'Interpreter', 'latex', 'fontsize',16)

figure;
plot(t, xo2, 'r');
xlabel('$t/s$', 'Interpreter', 'latex', 'fontsize',16)
ylabel('$x_o/m$', 'Interpreter', 'latex', 'fontsize',16)

result2 = [t; xf2; vf2; xo2; vo2]';
filename = '../../结果/result1-2.xlsx';
writematrix(result2,filename,'Sheet',1,'Range','A3:E900')

% 论文展现
p2_dis = deval(sol2,t_dis);
xf2_dis = p2_dis(1,:);
vf2_dis = p2_dis(2,:);
xo2_dis = p2_dis(3,:);
vo2_dis = p2_dis(4,:);
result2_dis = [t_dis; xf2_dis; vf2_dis; xo2_dis; vo2_dis]';
filename = '../../结果/result1-2_dis.xlsx';
writematrix(result2_dis,filename,'Sheet',1,'Range','A3:E7')
