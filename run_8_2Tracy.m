# Uebung 8 Aufgabe 2

function [] = run_8_2()
  
  n = 3;
  a = 0;
  b = 1;
  m = 1000;
  y_test1 = gauss_quad(@testfunction1, a, b, n, m)
  y_test1 = gauss_quad(@testfunction2, a, b, n, m)
  ## Tests geben nicht die richtigen Werte 
  #  --> Gewichte und Stuetzstellen checken!
end

function y = gauss_quad(fun, a, b, n, m)
# Berechnung eines Integrals mit Gauss-Legendre-Integration
# fun: Zu integrierende Funktion
# a  : Anfangswert des Integrals
# b  : Endwert des Integrals
# n  : Polynomgrad
# m  : Anzahl der Teilintervalle

## Stuetzstellen und Gewichte in Abh. von n (Wikipedia)
  switch n
    case 1
      mu = [2];
      knotes = [0];
    case 2
      mu = [1, 1];
      knotes = [-sqrt(1/3), sqrt(1/3)];
    case 3
      mu = [5/9, 8/9, 5/9];
      knotes = [-sqrt(3/5), 0, sqrt(3/5)];
    case 4
      mu = [ (18-sqrt(30))/36, (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36,]; 
      knotes =[ -sqrt( 3/7 + 2/7 * sqrt(6/5)), -sqrt( 3/7 - 2/7 * sqrt(6/5)), sqrt( 3/7 - 2/7 * sqrt(6/5)), sqrt( 3/7 + 2/7 * sqrt(6/5)),]; 
    otherwise
      disp('Polynomgrad soll kleiner als 5 sein!')
  end
  
  h = (b-a)/(m+1);
  for k = 1:m+1
    z(k) = a + k*h;
  end
  
  # Berechne x
  x_ik = zeros(n,m+1);
  
  for i = 1:n
    for k = 2:m+1
      x_ik(i,k) = z(k-1) + (h /(b-a)) *(knotes(i) - a);
    end
  end 
  
  # Berechne y
  y = 0;
  for k = 1:m+1
    for i = 1:n
      y = y +  mu(i) * fun(x_ik(i,k)) * h;
    end
  end  
  
end

## Test Funtkionen
function y = testfunction1(x)
  y = x;
end

function y = testfunction2(x)
  y = x**3;
end


