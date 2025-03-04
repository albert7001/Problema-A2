close all 
clc
clear

% Longitud de la red cuadrada
N = 100 ;

% Número de iteraciones
iter = 2000000 ; 


% Estado inicial del sistema (x es la red en representación matricial)
x = rand (N,N) - 0.5; % distribución uniforme entre -0.5 y 0.5
x = sign (x); % cambio del valor por su signo (matriz de 1 y -1 uniforme)

 
% Parámetros
 eps = 20 ; % Energía de acoplamiento
 h = 0;  % Campo magnético
 T = 0;  % Temperatura


% Ejecución de la función
[M] = monte_carlo(x,iter,N,h,eps,T);

M_rel = M/N^2 % magnetización relativa por unidad de espín



function [mag] = monte_carlo (x,iter,N,h,eps,T)

    figure
    imagesc(x) % matriz a imagen (1 blanco, -1 negro)
    colormap('gray')
    xlim ([1 N]) 
    ylim ([1 N])
    xlabel('x')
    ylabel('y')
    set(gcf,'Color','w')
    hold on 


    % Algoritmo de Metropolis
    for n = 1:iter 
  
        index = randi(N*N); % Se elige un nº al azar entre todos los espines 
        
        [i,j] = ind2sub([N N],index)  ; % Se traduce la posición lineal a la posición en la matriz
        
                                     
        % mod gives the report of the division between numerator and denominator
        % it's needed to take in account boundary condition
        
        dx = sub2ind([N N],i,(mod(j,N)+1)) ;  % índice del vecino derecho
        sx = sub2ind([N N],i,mod(j-2,N)+1) ;  % índice del vecino izquierdo
        up =sub2ind ([N N],mod(i-2,N)+1,j) ; % índice del vecino de arriba
        dn = sub2ind([N N],mod(i,N)+1,j) ;    % índice del vecino de abajo
        
        
        % Suma del valor del espín de los 4 próximos vecinos (los que interaccionan con el i,j)
        neigh = x(sx)+x(dx)+x(up)+x(dn)  ; 
        
        % Energía del espín (con signo cambiado)
        dE = eps*(x(index)*neigh) + h*(x(index)) ; 
        
        
        % Factor de Boltzmann
        prob = exp(-dE/T); 
        
        % Si la energía es positiva o por azar (controlado por el factor de Boltzmann) el espín cambia de sentido
        if dE <=0 || rand() <= prob
        
              x(index) = -x(index) ; 
    
        end
        
        f = 2000 ; % nº de iteraciones hasta mostrar la imagen
        
        if  sum(ismember(1:f:iter,n)) 
    
            imagesc(x)
            pause(0.0000000001) % tiempo de espera entre imágenes
    
       end   
 
   end 

    mag = sum(x(:)==1)-sum(x(:)==-1); % Cálculo de la magnetización final

end 
