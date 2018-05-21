%Numerik SoSe18 FU Berlin
%Aufgabenblatt 5
%David und Tracy

%Aufgabe 3 a) 
%Implementierung von AitkenNeville:

function [] = run_5_3()
  I = [-5,5];                     % Intervall
  n = 10;                         % Anzahl Stuetzstellen
  step = abs((I(2)-I(1))/n);
  x = zeros(1,n);
  x(1) = I(1);
  
  for i = 2:n
    x(i) = x(i-1) + step;         %Stuetzstellenwerte
  end
  
  % Funktionen zum Testen
  f1 = transpose(fOption_5_3(x,1));
  f2 = transpose(fOption_5_3(x,2));
  
  pointsAD_f1 = InterpolationAD(f1,I,n,x);
  pointsTscheby_f1 = InterpolationTscheby(f1,I,n);
  pointsAD_f2 = InterpolationAD(f2,I,n,x);
  pointsTscheby_f2 = InterpolationTscheby(f2,I,n);

  %{fplot(@(x) exp(-x), [x(1),x(end)]);
  hold on 
  plot(x, pointsAD_f1, '--b');
  plot(x, pointsTscheby_f1, '--r')
  hold off

  fplot(@(x) atan(x), [x(1),x(end)]);
  hold on 
  plot(x, pointsAD_f2, '--b');
  plot(x, pointsTscheby_f2, '--r')
  hold off
  %}
end

function v = AitkenNeville(x,y,u,n)
  % RUN_5_3 run Aufgabe 3  von Blatt 5
  %Initialisiere:
  AN = zeros(n,n);
  AN(:,1) = y(:);
  
  %Fuellen der Matrix nach dem AN Schema
  for i = 2:n
    j = 2;
    while (j <= n) && (j <= i) 
        AN(i,j)= (AN(i,j-1)-AN(i-1,j-1))/(x(i)-x(i-j+1));
        j = j+1;
    end
  end
  
  %Initialisierung zur Berechnung der Newtonschen Polynome:
  coeff = zeros(1,n);
  k = 0;
  for i = 1:n
    coeff(i)=AN(i,i); 
    v = coeff(i) * (u**k);
    k = k + 1;
  end
    
end 


%Aufgabe 3 b)
%Interpolation fuer Aequidistante Stuetzstellen:

function f_interpolAD = InterpolationAD(f,I,n,x)

  f_interpolAD = zeros(1:n,1);
  
  for i = 1:n
    f_interpolAD(i) = AitkenNeville(x,f,i,n);
  end
end

%Aufgabe 3c)
function f_interpolTscheby = InterpolationTscheby(f,I,n)
  x = zeros(1:n,1);
  for i = 1:n
     x(i)=cos((2*(n-i)+1)*pi/(2*n+2));
  end

  f_interpolTscheby = zeros(1:n,1);

  for i = 1:n
    f_interpolTscheby(i) = AitkenNeville(x,f,i,n);
  end
end



