%Numerik SoSe18 FU Berlin
%Aufgabenblatt 5
%David und Tracy

%Aufgabe 3 a) 
%Implementierung von AitkenNeville:

function [] = run_5_3()
I = [-5,5]; 

InterpolationAD(f1,I,5)
InterpolationTscheby(f1,I,5)
InterpolationAD(f2,I,5)
InterpolationTscheby(f2,I,5)
end

function v = AitkenNeville(x,y,u)
  % RUN_5_3 run Aufgabe 3  von Blatt 5
  %Initialisiere:
  n = length(x);
  AN = zeros(n);
  AN(1,:) = y;
  
  %Füllen der Matrix nach dem AN Schema
  for i = 2:n
    for j = 1:n
      if i <= j
        AN(i,j)= (AN(i-1,j)-AN(i-1,j-1))/(x(j)-x(j-i+1));
      end
    end
  end
  
  %Initialisierung zur berechnung der Newtonschen Polynome:
  coeff = [1:n];
  polynom = coeff;
  
  
  for k = 1:n
    coeff(k)=AN(k,k);
    polynom(k+1) = polynom(k)*(u-x(k)); 
  end
  v = dot(coeff,polynom);
end

%Aufgabe 3 b)
%Interpolation für äquidistante Stützstellen:

function [] = InterpolationAD(f,I,n)
  x = [I(1):(I(2)-I(1))/n:I(2)];
  y = f(x);
  f_interpol = (1:n+1);
  
  for l = 1:n+1
    f_interpol(l) = AitkenNeville(x,y,x(l));
  end
  
  plot(x,f_interpol)
end

%Aufgabe 3c)
function [] = InterpolationTscheby(f,I,n)
  x=(0:1:n);
  for h = 1:n+1
     x(h)=cos((2*(h-1)+1)*pi/(2*n+2));
end

  y = f(x);
  f_interpol = (1:n+1);

  for l = 1:n+1
    f_interpol(l) = AitkenNeville(x,y,x(l));
  end
  
  plot(x,f_interpol)
end

%Aufgabe 3e):
%Test der Funktionen:

function y = f1(x)
 y = exp(-x);
end

function y = f2(x)
y = atan(x);
end


