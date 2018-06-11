%Numerik I SoSe 18 
%David und Tracy
%Aufgabe 2


function [Integral] = run_8_2()
%Test der Gauss Quadratur am Polynom x^3 im Intervall [0,1] 
[Integral] = gauss_quad(@func,0,1,3,20);

%Anwendung der Funktion gauß_quadratur zu f(x) = 2*log(x/2+1)/(x^2+4) in [0,2]
%Dabei wird die Quadratur für verschiedene Werte n =1,2,3,4 und m = 1,2,...,10 durchgeführt.
%Anschließend wird der Fehler logarithmisch über m Aufgetragen und für alle n geplottet.
%Daraus lässt sich die Konvergenzordnung ablesen.(Aufgabe 2b)

%Initialisieren
I = zeros(1,10);
m = 1:10;

%For Schleifen zur variation von n, m :
for s = 1:4;
for i = 1:10
  I(i)=gauss_quad(@f,0,2,s,i);
end
%bestimmung und plot des Fehlers:
Ierr = abs(I-log(2)*pi/8);
figure
loglog(m,Ierr)
xlabel('ln(m)')
ylabel('ln(Fehler)')
title(strcat('Logarithmische Auftrageung der Fehler zu n =', num2str(s)))
end

end
%Implementierung der Gaußquuadratur mit Koeffizienten aus Bronstein:(Aufgabe 2a)
function [y] = gauss_quad(f,a,b,n,m)
x = linspace(a,b,m+1);
xi= zeros(4);
ai = zeros(4);

xi(1,1)=0;
xi(1:2,2) = [-sqrt(1/3) sqrt(1/3)];
xi(1:3,3) = [-sqrt(3/5) 0 sqrt(3/5)];
xi(1:4,4) = [-sqrt(3/7+(2/7)*sqrt(6/5)) -sqrt(3/7-(2/7)*sqrt(6/5)) sqrt(3/7-(2/7)*sqrt(6/5)) sqrt(3/7+(2/7)*sqrt(6/5))];

ai(1,1) = 2;
ai(1:2,2)=[1 1];
ai(1:3,3)=[5/9 8/9 5/9];
ai(1:4,4)=[(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36];


int = 1:m;
for i = 1:m
    yi = ((x(i+1)-x(i))/2)*xi+(x(i+1)+x(i))/2;
    aai = ((x(i+1)-x(i))/2)*ai;
    int(i)= dot(f(yi(1:n,n)),aai(1:n,n));
end
y = sum(int);
end 

function f = func(x)
f = x.^3;
end

function L = f(x)
  L = 2*log(x/2+1)./(x.^2+4);
end


