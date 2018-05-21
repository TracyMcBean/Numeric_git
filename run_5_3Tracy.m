% Uebung 5, Aufgabe 3
% Tracy, David
% Abgabe 22.05.2018

function c = run_5_3()
% RUN_5_3 run task 3 from homework 5.
  n = 10;
  step = 10/n;
  x = [-5:step:5];

  f_1 = exp(-x);
  f_2 = atan(x);
  
  % Test 1
    fprintf(" \n Test 1 with function 1: exp(-x) \n")

  for i = -5:1:5
    c = AitkenNeville(x,f_1,i);
    disp(c);
  end
  % Test 2
  fprintf(" \n Test 2 with function 2: arctan(x) \n")
  for i = -5:1:5
    c = AitkenNeville(x,f_2,i);
    disp(c);
  end 
   

end

function v = AitkenNeville(x,y,u) 
% AITKENNEVILLE run the Aitken Neville scheme. 
%   V = AitkenNeville(x,y,u) calculates function value for an interpolated 
%       polynome p(x)=y at position u with given grid points x.
  l = length(x);
  % Initialize cooeficient vector for polynome
  p = zeros(1,l);
  % Initialize matrix for Neville Table
  f = zeros(l,l);
  % Aitken Neville scheme
  
  f(1,1) = y(1);
  f(:,1) = y(:);
  j   = 2; 
  for i = 2:l
    while (j <= i) && (j <= l) 
      f(i,j) = (f(i,j-1)-f(i-1,j-1))/(f(i,1) - f(i-j+1,1));
      j = j + 1;
    end
  end 
  
  for i = 1:l
    p(i) = f(i,i);
  end 
  
  % for loop to calculate the result
  % k as exponential 
  k = 0;
  for i = 1:length(x)
    v = p(i) * (u**k);
    k = k + 1;
  end 
  
end