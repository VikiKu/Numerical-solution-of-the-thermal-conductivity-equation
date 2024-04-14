close all
clear all
%% объявление переменных
%шаг
dx = 0.1;
dt = 0.005; % а)0.01 б)0.005
%сетка
x = 0:dx:10; 
t = 0:dt:1;
%задаю массив температуры
T = zeros(length(x), length(t));
%задаю граничные условия
T(1, :) = 0;
T(end, :) = 0;
%задаю начальные условия
T0 = 2.27;
x0 = 5;
T(:, 1) = T0*(x - x0).^2.*exp(-(x - x0).^2);

% шесть моментов времени, для которых просят построить графики Т(х)
i1 = find(t == 0);
i2 = find(t == 0.1);
i3 = find(t == 0.2);
i4 = find(t == 0.3);
i5 = find(t == 0.5);
i6 = find(t == 1);
I = [i1, i2, i3, i4, i5, i6];

%% явная схема %x -> j, t -> n
alpha = dt/dx^2;
for it = 1:length(t) - 1
 for ix = 2:length(x) - 1
 T(ix, it + 1) = T(ix, it) + alpha*(T(ix - 1 , it)...
 - 2*T(ix, it) + T(ix + 1, it)); %шаблон явной схемы
 end
end
%% построение графиков явной схемы решения
figure(1)
plot(x, T(:, I), 'LineWidth', 2)
set(gca,'FontSize',20,'FontName','Times New Roman')
legend('t = 0', 't = 0.1','t = 0.2', 't = 0.3', 't = 0.5', 't = 1')
xlabel('x')
ylabel('T(x, t)')
title('\Delta x = 0.1, \Delta t = 0.01') %тут исправлять при расчете А и Б
%%
figure(2)
for k = 1:6
subplot(2, 3, k)
plot(x, T(:, I(k)));
title(['t = ', num2str(t(I(k)))]);
xlabel('x')
ylabel('T(x, t)')
end
%% схема кранка-николсона
% x -> j, t -> n
a = 0*x; 
a(end) = 0;
b = 0*x;
b(end) = 0;
f = 0*x;
f(end) = 0;
% граничные условия
T2 = 0*T;
T2(1, :) = 0;
T2(end, :) = 0;
% начальные условия
T2(:, 1) = T0*(x - x0).^2.*exp(-(x - x0).^2);

%%
for n = 1:length(t) - 1
 % расчет коэффициентов решения
 for j = length(x)-1:-1:2
 a(j - 1) = 1/(2*(1 + 1/alpha) - a(j));
 f(j) = T2(j + 1, n) + 2*(1/alpha -1)*T2(j, n) + T2(j - 1, n);
 b(j - 1) = (f(j) + b(j))/(2*(1/alpha + 1) - a(j));
 end
 %решение
 for j = 1:length(x) - 1
 T2(j + 1, n + 1) = a(j)*T2(j, n + 1) + b(j);
 end
end

%% построение графиков FN
figure(3)
plot( x, T2(:, I), 'LineWidth', 2)
legend('t = 0', 't = 0.1','t = 0.2','t = 0.3', 't = 0.5', 't = 1')
xlabel('x')
ylabel('T(x, t)')
title('\Delta x = 0.1, \Delta t = 0.005')

%% сеточная диффузия
%шаг по времени
dt = 0.01;

f = (0:0.001:1)*pi/2; %k/k_n*pi/2

G1 = -log(abs(1 - 4*dt/dx/dx.*(sin(f)).^2))/pi; 
G2 = -log(abs((1 - 2*dt/dx/dx.*(sin(f)).^2)./((1 + 2*dt/dx/dx.*(sin(f)).^2))))/pi;
G3 = dt/pi* (2/dx.*f).^2; %(f.*2/dx).^2/pi; %dt/2* (pi^2/2/dx.*f).^2; 
%построение графиков
figure(4)
plot(f*2/pi, G1,f*2/pi, G2, f*2/pi, G3, 'LineWidth', 1.5);
legend('явная схема', 'схема Кранка-Николсона', 'исходное уравнение')
xlabel('k/k_N')
ylabel('\Gamma\Deltat/\pi')
ylim([-1,3])
title('\Delta t = 0.005')



