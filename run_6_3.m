%Aufgabenblatt 6 Aufgabe 3
%David und Tracy

%3 a)
%Test der unten definierten Funktion Plotspline an der Funktion func und
%Gitter.x
%Die Implementierung ist zudem Teil der Aufgabe

x = linspace(-1,1,5);
f = func(x);

erro = plotspline(x,@func);

%3 b)
%Darstellung der Fehler fuer verschiedene Gitter:
%Initialisieren:
n(1:10) = 2.^(1:10);
error = 1:10;
ddf_max = 1:10;

for i = 1:10
    x = linspace(-1,1,n(i)+1);
    error(i)= plotspline(x,@func);
    %Werte der zweiten Ableitung zur Fehlerabschaetzung:
    ddf = -50./(1+25*x.^2).^2 + (50*100.*x.^2)./(1+25*x.^2);
    ddf_max(i)= max(ddf);
end

lin(1:10) = (max(ddf_max).*(2./n).^2)./8;
figure
loglog(n,error)

%Plot des maximalen Fehlers nach Vorlesungsmitschrift: 
%err <= (||f''||*h^2)/8
hold on
loglog(n,lin)
hold off


function f = func(x)
f(:) = 1./(1+25*x.^2);
end

function [err] = plotspline(x,f)
%PLOTSPLINE Die Funktion plottet die Funktion f(x), sowie den
%interpolierten linearen Spline und gebit den Fehler in Maximumsnorm aus.
% Eingabe Parameter sind dabei die Funktion, sowie das gewaehlte Gitter.

%Initialisiere
Gitter = x;
f_x = f(x);

svec = zeros(1,(length(x)-1)*20);

fvec = zeros(1,(length(x)-1)*20);

y_ges = zeros(1,(length(x)-1)*20);

error = zeros(1,(length(x)-1)*20);

%For Schleife zum Auswerten der Funktion, beziehungsweise der Splines
for i = 1:length(x)-1
    
    y = linspace(Gitter(i),Gitter(i+1),20);
    y_ges(((i-1)*20+1):i*20) = y;
    m = (f_x(i+1)-f_x(i))/(Gitter(i+1)-Gitter(i));
    c = f_x(i)-m*Gitter(i);
    svec(((i-1)*20+1):i*20) = m*y+c ;
    fvec(((i-1)*20+1):i*20) = f(y);
    error = max(abs(svec(((i-1)*20+1):i*20)-fvec(((i-1)*20+1):i*20)));
   
end

%Darstellung der Funktion, sowie der Splines:
 plot(y_ges,svec)
 hold on
 plot(y_ges,fvec)
 hold off 
 
%Berechnung des maximalen errors.
err = max(error);
end



