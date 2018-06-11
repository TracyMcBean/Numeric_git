# Numerik 1, SS18
# Uebung 8 Aufgabe 1
# David und Tracy

function [] = run_8_1()
  
  printf("Aufgabe 1a:\n")
  # Trapezregel
  # (Test 1 fuer f(x)=x da Exaktheit fuer grad = 1 gilt)
  a = 0; 
  b = 1;
  n = 20;
  printf("Test mit y = x fuer Trapezregel \n")
  y_trapez1 = trapez(@testfunction1, a, b, n)
  printf("Test mit y = x^2 fuer Trapezregel \n")
  y_trapez2 = trapez(@testfunction2, a, b, n)
  
  ## Simpsonregel
  printf("Test mit y = x fuer Simpsonregel \n")
  y_simpson1 = simpson(@testfunction1, a, b, n)
  printf("Test mit y = x^2 fuer Simpsonregel \n")
  y_simpson2 = simpson(@testfunction2, a, b, n)
  
  # Aufgabe 1b
  ## Verwende exp(x) und sqrt
  a = -2; b = 2;
  i = 1;
  for n = 2:1:100
    steps(i) = n;
    # Test mit y = exp(x) fuer Trapezregel
    y_trapez_exp(i) = trapez(@expfunction, a, b, n);
    # Test mit y = sqrt(abs(x)) fuer Trapezregel 
    y_trapez_sqrt(i) = trapez(@rootfunction, a, b, n);
    i = i + 1;
  end
  
  error_trapez_exp = abs(y_trapez_exp(:) - (exp(2) - exp(-2)) );
  error_trapez_sqrt = abs(y_trapez_sqrt(:) - (2/3*(2)**3/2 - 2/3*(-2)**3/2) ); 
   
  a = -9; b = 9;
  i = 1;
  for n = 2:1:100
    #Test mit y = exp(x) fuer Simpsonregel 
    y_simpson_exp(i) = simpson(@expfunction, a, b, n);
    #Test mit y =sqrt(abs(x)) fuer Simpsonregel
    y_simpson_sqrt(i) = simpson(@rootfunction, a, b, n);
    i = i + 1;
  end
  
  error_simpson_exp = abs(y_simpson_exp(:) - (exp(9) - exp(-9)) );
  error_simpson_sqrt = abs(y_simpson_sqrt(:) - (2/3*(9)**3/2 - 2/3*(-9)**3/2) ); 

  ## Plot fuer Fehlerdarstellung
  figure
  hold on
  plot(steps(:), error_trapez_exp(:), 'LineWidth',2)
  plot(steps(:), error_trapez_sqrt(:), '--g','LineWidth',2)
  hold off
  
  figure
  plot(steps(:), error_simpson_exp(:), 'LineWidth',2)
  figure
  plot(steps(:), error_simpson_sqrt(:), '--g', 'LineWidth',2)
end

function y = trapez(fun, a, b, n)
# Berechnung eines Integrals anhand der Trapezregel
# fun: Zu integrierende Funktion
# a  : Anfangswert des Integrals
# b  : Endwert des Integrals
# n  : Anzahl der Teilintervalle
  if a > b
    printf("Warnung: Anfangswert muss kleiner als Endwert sein!")  
  elseif n < 1 
    printf("Warnung: n muss groeßer als 0 sein!")
  end
  
  # Teilintervalle definieren
  I_x = linspace(a,b,n);
  # Berechne die Summe der Integralle über die Teilintervalle mit Trapezregel
  y = 0;
  for i = 1 : n-1
     y = y +( ( I_x(i+1) - I_x(i) ) * (fun(I_x(i)) + fun(I_x(i+1)) ) * 1/2 );
  end
  
end

function y = simpson(fun,a,b,n)
# Berechnung eines Integrals anhand der Simpsonregel
# fun: Zu integrierende Funktion
# a  : Anfangswert des Integrals
# b  : Endwert des Integrals
# n  : Anzahl der Teilintervalle
  if a > b
    printf("Warnung: Anfangswert muss kleiner als Endwert sein!")  
  elseif n < 1 
    printf("Warnung: n muss groeßer als 0 sein!")
  end
  
  # Teilintervalle definieren
  I_x = linspace(a,b,n);
  # Berechne die Teilintervalle
  y = 0;
  for i = 1 : n-1
    c = fun( (I_x(i)+ I_x(i+1))/2 );
    y = y + ( I_x(i+1) - I_x(i) )/6 * ( fun(I_x(i)) + 4*c + fun(I_x(i+1)) );
  end
  
end

## Test Funtkionen
function y = testfunction1(x)
  y = x;
end

function y = testfunction2(x)
  y = x**2;
end


### Funktione fuer b
function y = expfunction(x)
  y = exp(x);
end

function y = rootfunction(x)
  y = sqrt(abs(x));
end
