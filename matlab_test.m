%Teste 1

F = @(t, x) 1 + (x - t)^2;
ts = 1.05:0.01:3;
figure(1)
ode45(F, ts, -18.95);

%Teste 2

Sistema = @(t,x) [ (-2*x(1) -1*x(2) -1*x(3) -2*x(4));
    (x(1) -2*x(2) +2*x(3) -1*x(4));
    (-1*x(1) -2*x(2) -2*x(3) -1*x(4));
    (2*x(1) -1*x(2) +1*x(3) -2*x(4))
    ];
t_span = 0:0.01:2;
[t,x] = ode45(Sistema,t_span,[1;1;1;-1])
figure(2)
plot(t,x)

