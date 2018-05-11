% Uebung 4, Numerik 1
% David, Tracy
% Aufgabe 1

function result = run_1_4()
    %RUN_1_4() loest lineares Ausgleichsproblem fuer verschieden Gitter
    
    disp('Test 1, Gitter mit m = 2 und n = 1:')
    result = givens_rotation(2,1);
    fprintf('Test result: %s \n', result)
    
    disp('----------------------------------------------')
    disp('Test 2, Gitter mit m = 4 und n = 2:')
    result =  givens_rotation(4,2);
    fprintf('Test result:\n')
    disp(result)
    
    disp('----------------------------------------------')
    disp('Test 3, Gitter mit m = 10 und n = 3:')
    result =  givens_rotation(10,3);
    fprintf('Test result:\n')
    disp(result)
    
end %function

function result = givens_rotation(m, n)
    %GIVENS_ROTATION() fuehrt Givens Rotation aus um lin. Ausgleichsproblem
    %zu loesen.
    % n: Grad vom Polynom
    % m: Gitteraufloesung
    
    % Initialisiere
    x = zeros(m, 1);
    b = zeros(m, 1);
    xcooef = zeros(m,1);  % fuer matrix A gebraucht
    A = zeros(m,n+1);
    x(m)= 2*pi;
    xcooef(:) = 1;
    
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
    for i = 1:(n+1)
        A(:, i) = xcooef;
        xcooef = xcooef.* x;
    end 
    
    % Erstelle 3D-Matrix mit Einheitsmatrix
    Q_single = eye(m);
    % TO DO number of Qs
    for k = 1:(m*n)
        Q(:,:,k) = Q_single;
    end 
    
    % Berechne die Transformationsmatrix
    count = 1;
    for i = 1:(n+1)
        a = A(:,i);
        l = 1;
        k = 1;
        while (l< m)
            a_l = a(l);
            a_k = a(l+1);
            r = sqrt(a_l^2 + a_k^2);
            s = a_k / r;
            c = a_l / r;        
            
            Q(l,k,count) = -s;
            Q(l,l,count) = c;
            Q(k,l,count) = s;
            Q(k,k,count) = c;
            
            count = count +1;
            l = l + 1;
        end   
    end
  
    Qnew = Q(:,:,1);
    for i = 2:(count-1)
        Qnew = Qnew * Q(:,:,i);
    end
    QT = transpose(Qnew);
    
    bEnde = QT * b;
    b_1 = bEnde(1:n+1);
    
    result = b_1;
    
end %function