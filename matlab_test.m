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
[t,x] = ode45(Sistema,t_span,[1;1;1;-1]);
figure(2);
plot(t,x);

%Teste 3

n = 5;
A = zeros(n,n);

for i = 1:n
    for j = 1:n
        if i == j
            A(i,j) = -2;
        elseif j == i + 1
                A(i,j) = 1;
        else
            A(i,j) = 0;
        end
    end
end

A;

X0 = zeros(n,1);

for i = 1:n
    y = i/(n+1);
    X0(i) = sin(pi*y) + sin(n*pi*y);
end

X0;

 % lembrando n = 5
sistema3 = @(t,x) [(A(1,1)*x(1)+A(1,2)*x(2)+A(1,3)*x(3)+A(1,4)*x(4)+A(1,5)*x(5));
    (A(2,1)*x(1)+A(2,2)*x(2)+A(2,3)*x(3)+A(2,4)*x(4)+A(2,5)*x(5));
    (A(3,1)*x(1)+A(3,2)*x(2)+A(3,3)*x(3)+A(3,4)*x(4)+A(3,5)*x(5));
    (A(4,1)*x(1)+A(4,2)*x(2)+A(4,3)*x(3)+A(4,4)*x(4)+A(4,5)*x(5));
    (A(5,1)*x(1)+A(5,2)*x(2)+A(5,3)*x(3)+A(5,4)*x(4)+A(5,5)*x(5));
    ];
tspan = 0:0.01:2;
[t,x] = ode45(sistema3,tspan,X0);
figure(3)
plot(t,x);

% P�ndulo
g = 9.8;
l = 0.8;
w = sqrt(g/l);
u = 2*w - 0.5;
sistema = @(t,x) [(x(2)),
                  (-(w^2)*sin(x(1)) - u*x(2))]
tspan = 0.0:0.01:5;
[t,x] = ode45(sistema,tspan,[(5*pi)/2, 3]);
figure(4)
plot(t,x)

