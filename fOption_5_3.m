% Datei mit Funktionen
function f = fOption_5_3(x, f_opt)
  % Falls keine Funktionsnummer angeben gilt default = exp(-x)
  if nargin < 2
    f_opt = 1;
  end 

  if f_opt == 1
    f =  exp(-x);
  elseif f_opt == 2
    f = atan(x);
  else 
    fprintf("Nur 2 Beispiel Funktionen vorhanden.")
  end
end