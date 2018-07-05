%Tracy und David
%Zettel 12

%Aufgabe 1)

tau = 0.01;
gitter = 0:tau:1000;
A = [0 -1 ; 1 0];
x0 = [0 1]';

ex = expl(A,gitter,x0);
im = impl(A,gitter,x0);

scatter(ex(1,:), ex(2,:))
hold on
scatter(im(1,:),im(2,:))
hold off

%Aufgabe 2)
A = [0 0 0 0 ; 0.5 0 0 0 ; 0 0.5 0 0 ; 0 0 0.5 0];
b = [1/6 1/3 1/3 1/6]';
tspan = [0 2];
x0 = -pi/4;
fex = @(t)atan(tan(x0)+100*t);

%a)
f = @(x)100/(1+(tan(x))^2);
N = [2^4 2^5 2^6 2^7 2^8 2^9];
[t1,x1,c1] = RungeKuttaEx(f,tspan,x0,N(1),b,A);
[t2,x2,c2] = RungeKuttaEx(f,tspan,x0,N(2),b,A);
[t3,x3,c3] = RungeKuttaEx(f,tspan,x0,N(3),b,A);
[t4,x4,c4] = RungeKuttaEx(f,tspan,x0,N(4),b,A);
[t5,x5,c5] = RungeKuttaEx(f,tspan,x0,N(5),b,A);
[t6,x6,c6] = RungeKuttaEx(f,tspan,x0,N(6),b,A);

figure
plot(t6,fex(t6))
hold on
plot(t1,x1)
plot(t2,x2)
plot(t3,x3)
plot(t4,x4)
plot(t5,x5)
plot(t6,x6)
hold off
legend('exakt','N = 2^4','N = 2^5','N = 2^6','N = 2^7','N = 2^8','N = 2^9')

%b)

gamma = [1 10 100 1000];
e = zeros(6,4);
cges = zeros(6,4);

%Erste Schleife ueber die Gamma-Werte
for l = 1:4
   fun = @(x)gamma(l)/(1+(tan(x))^2);
   fex = @(t)atan(tan(x0)+gamma(l)*t);
   %Zweite Schleife ueber die N-Werte
   for m = 1:6
       [t,x,c] = RungeKuttaEx(fun,tspan,x0,N(m),b,A);
       e(m,l)= max(abs(fex(t)-x));
       cges(m,l)=c;
   end
%LogLog-Plot der Fehler ueber den f-Auswertungen
figure
loglog(cges(:,l),e(:,l))
title(['Gamma =' num2str(gamma(l))])
xlabel('f-Auswertungen')
ylabel('Fehler')
end

%Implementierung explizites Runge-Kutta-Verfahren:

function [t,x,c] = RungeKuttaEx(f,tspan,x0,N,b,A)
s = length(b);
k=b;
x = 1:(N+1);
tau = (tspan(2)-tspan(1))/N;
x(1) = x0;
c = 0;
 for j = 2:(N+1)
    k(1)=f(x(j-1));
    c=c+1;
    for i = 2:s
        k(i)= f(x(j-1)+tau*dot(A(i,1:(i-1))',k(1:(i-1))))';
        c=c+1;
    end
    x(j) = x(j-1)+tau*dot(b,k');
 end
 t = tspan(1):tau:tspan(2);
end


%Implementierung des expliziten Eulers
function y = expl(A,gitter,x0)
  y = zeros(length(x0),length(gitter));
  y(:,1) = x0;
  for i = 2:length(gitter) 
      y(:,i)= y(:,i-1)+(gitter(i)-gitter(i-1))*A*y(:,i-1);
  end
end

%Implementierung des impliziten Eulers
function y = impl(A,gitter,x0)
  y = zeros(length(x0),length(gitter));
  y(:,1) = x0;
  for i = 2:length(gitter)
      y(:,i)= (eye(length(x0))-A)\y(:,i-1);
  end
end

