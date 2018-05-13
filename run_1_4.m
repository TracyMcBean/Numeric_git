% Uebung 4, Numerik 1
% David, Tracy
% Aufgabe 1

function result = run_1_4()
    %RUN_1_4() loest lineares Ausgleichsproblem fuer verschieden Gitter
    
    disp('Test 1, Gitter mit m =3 und Polynom Grad n = 1:')
    result = givens_rotation(3,1);
    fprintf('Ergebnis fuer u:\n')
    disp(result)
    
    disp('----------------------------------------------')
    disp('Test 2, Gitter mit m = 10 und Polynom Grad n = 3:')
    result =  givens_rotation(10,3);
    fprintf('Ergebnis fuer u:\n')
    disp(result)
    
    disp('----------------------------------------------')
    disp('Test 3, Gitter mit m = 4 und Polynom Grad n = 4:')
    result =  givens_rotation(4,4);
    fprintf('Ergebnis fuer u:\n')
    disp(result)
    
end %function

function result = givens_rotation(m, n)
    %GIVENS_ROTATION() fuehrt Givens Rotation aus um lin. Ausgleichsproblem
    %zu loesen.
    % n: Grad vom Polynom
    % m: Gitteraufloesung
    
    % Fange Fehler ab wenn m zu niedrig ist
    if (m < 3)
      disp('Error: m must be bigger than 2!')
      exit;
    elseif (m < n)
      disp('Error: m must be bigger than n!')
      exit;
    end
    
    % Initialisiere
    x = zeros(m, 1);
    b = zeros(m, 1);
    xcooef = zeros(m,1);  % fuer matrix A gebraucht
    A = zeros(m,n);
    x(m)= 2*pi;
    xcooef(:) = 1;
    n = n+1;    % Weil polynom vom Grad n zusaetzlich x^0 beinhaltet
    
    % Falls mehr als 2 Gitterpunkte setze die Abstaende entsprechend.
    if ( m > 2)
        step = 2*pi / m;
        for i = 2:(m-1)
            x(i) = x(i-1) + step;
        end %for
    end % if
 
    for i = 1:m
        b(i) = sin(x(i));
    end 
    
    % Matrix A auffuellen mit t
    for i = 1:(n)
        A(:, i) = xcooef;
        xcooef = xcooef.* x;
    end 
    
    % Erstelle 3D-Matrix mit Einheitsmatrix die spaeter aufgefuellt wird.
    Q_single = eye(m);
    for k = 1:(m*n-m)
        Q(:,:,k) = Q_single;
    end 
    
    % Berechne die Transformationsmatrix
    count = 1;
    for i = 1:(n)
        a = A(:,i);
        l = i;
        k = l+1;
        while (l< m)
            a_l = a(l);
            a_k = a(k);
            r = sqrt(a_l^2 + a_k^2);
            s = a_k / r;
            c = a_l / r;        
            
            Q(l,k,count) = s;
            Q(l,l,count) = c;
            Q(k,l,count) = -s;
            Q(k,k,count) = c;
            
            count = count +1;
            l = l + 1;
            k = k + 1;
        end   
    end
    
    % Multipliziere Transformationsmatrizen miteinander.
    Qnew = Q(:,:, count-1 );
    for i = count-1 : -1 : 1
        Qnew = Qnew * Q(:,:,i);
    end

    QT = transpose(Qnew);
    
    % Berechne neuen b Vector
    bEnde = QT * b;
    b_1 = bEnde(1:n);
    
    
    QTtimesA = QT * A;
    
    % Berechne R 
    R = QTtimesA(1:n, 1:n);
    
    % Berechne x um damit u zu bestimmen
    x = inv(R) * b_1;
    fprintf('Berechnetes x: \n')
    disp(x);
    
    u = A*x;
    result = u;
    
end %function