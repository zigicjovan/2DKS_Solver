function result = multiply2D(u,v,space)

% Multiply two 2D vector u,v using DEALIASING
% 'space' refers to the space on which the vectors are defined:
% space = 'fourier' or space = 'real'
    
    Nx = size(u,1);
    Ny = size(u,2);
    Nx2 = floor(Nx/2);
    Ny2 = floor(Ny/2);
    
    if mod(Nx2,2) ~= 0
        Nx2 = Nx2 - 1;
    end
    if mod(Ny2,2) ~= 0
        Ny2 = Ny2 - 1;
    end
    
    Pad = ones(Nx,Ny);
    myfilter = [ Pad( 1:Ny2/2 , 1:Nx ) ;
                 Pad( 1:Ny2 , 1:Nx2/2 ), zeros(Ny2,Nx-Nx2), Pad( 1:Ny2 , 1:Nx2/2 ) ;
                 Pad( 1:Ny2/2 , 1:Nx )]; 
    
    switch space
        case 'real'
            u_hat = myfilter.*fft2(u);
            v_hat = myfilter.*fft2(v);
            
            u = real(ifft2(u_hat));
            v = real(ifft2(v_hat));
            
            aux = myfilter.*fft2(u.*v);
            
            result = real(ifft2(aux));            
        case 'fourier'
            u_hat = myfilter.*u;
            v_hat = myfilter.*v;
            
            u = real(ifft2(u_hat));
            v = real(ifft2(v_hat));
            result = myfilter.*fft2(u.*v);
    end

return