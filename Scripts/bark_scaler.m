f = 1:500:8e3;

[b1, c1] = v_frq2bark(f, 'lh');
[b2, c2] = v_frq2bark(f, 'LH');

v_frq2mel(f);

% figure();
% stem(b1);
% hold on;
% stem(b2);
% grid on;
% legend on;

% figure();
% stem(c);
% grid on;
% legend on;

a = 10000;

x = 20:20e3;
t1 = log(x);
t2 = a*x.^(1/a)-a;
% t3 = x.^2 / 2;

t4 = x.^(1/a);

figure();
plot(t1);
hold on;
grid on;
plot(t2);
% plot(t3);
% plot(t4);
legend on;