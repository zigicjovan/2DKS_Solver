function result = multiply2D(u,v,space)

    switch space
        case 'fourier'
    
            % temporarily here to avoid dealiasing

            % Multiply two 2D vector u,v using DEALIASING
            % 'space' refers to the space on which the vectors are defined:
            % space = 'fourier' or space = 'real'
            
            Nx = size(u,1);
            Ny = size(u,2);
            
            % 1/2 dealiasing
            %{
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
            %
            % 2/3 dealiasing
            Nx3 = floor(Nx/3*2);
            Ny3 = floor(Ny/3*2);
        
            Pad = ones(Nx,Ny);
            myfilter = [ Pad( 1:ceil(Ny3/2)-1 , 1:Nx ) ;
                         Pad( 1:ceil(Ny3/2)+1 , 1:ceil(Nx3/2) ), zeros(ceil(Ny3/2)+1,Nx-Nx3-1), Pad( 1:ceil(Ny3/2)+1 , 1:ceil(Ny3/2) ) ;
                         Pad( 1:ceil(Ny3/2)-1 , 1:Nx )]; 
            %}
            % Gaussian spectral filter
            %{
            k1_n_pts = [ 0 : Nx/2-1 , 0 , -Nx/2+1 : -1]; 
            k2_n_pts = [ 0 : Ny/2-1 , 0 , -Ny/2+1 : -1]; 
            [ k1_n , k2_n ] = meshgrid(k1_n_pts,k2_n_pts); % 2-dimensional grid
            %k_n = abs(k1_n) + abs(k2_n);
            k_n = sqrt(abs(k1_n).^2 + abs(k2_n).^2);
            %}
            k1_l_pts = [0:Nx/2 -Nx/2+1:-1];
            k2_l_pts = [0:Ny/2 -Ny/2+1:-1];
            [ k1_l , k2_l ] = meshgrid(k1_l_pts,k2_l_pts); % 2-dimensional grid
            %k_n = abs(k1_l) + abs(k2_l);
            k_n = sqrt(abs(k1_l).^2 + abs(k2_l).^2);
            k_cut = sqrt((512/2 + 1).^2 + (512/2 + 1).^2)*(2/3);
            %k_cut = sqrt((512/2 + 1).^2 + (512/2 + 1).^2)*(2/3);
            %
        
            myfilter = NaN(Nx,Ny);
            for i = 1:Nx
                for j = 1:Ny
                    %myfilter(i,j) = exp(-36*(k_n(i,j)/(Nx/3*2)).^36);
                    myfilter(i,j) = exp(-36*(k_n(i,j)/k_cut).^36);
                end
            end
            %{
            % No dealiasing
            Pad = ones(Nx,Ny); 
            myfilter = Pad;
            %}

            u_hat = u;
            v_hat = v;
            
            u = real(ifft2(u_hat));
            v = real(ifft2(v_hat));
            result = myfilter.*fft2(u.*v);
            result(abs(result) < 1e-16) = 0;


        case 'real'
            u_hat = myfilter.*fft2(u);
            v_hat = myfilter.*fft2(v);
            
            u = real(ifft2(u_hat));
            v = real(ifft2(v_hat));
            
            aux = myfilter.*fft2(u.*v);
            
            result = real(ifft2(aux));            
        case 'dealias'

            Nx = size(u,1);
            Ny = size(u,2);

            k1_l_pts = [0:Nx/2 -Nx/2+1:-1];
            k2_l_pts = [0:Ny/2 -Ny/2+1:-1];
            [ k1_l , k2_l ] = meshgrid(k1_l_pts,k2_l_pts); % 2-dimensional grid
            %k_n = abs(k1_l) + abs(k2_l);
            k_n = sqrt(abs(k1_l).^2 + abs(k2_l).^2);
            k_cut = sqrt((512/2 + 1).^2 + (512/2 + 1).^2)*(2/3);

            myfilter = exp(-36*(k_n/k_cut).^36);

            result = myfilter.*u; 
            result(abs(result) < 1e-16) = 0;
        case 'real2'
            u_hat = myfilter.*fft2(u);
            v_hat = myfilter.*fft2(v);
            
            u = real(ifft2(u_hat));
            v = real(ifft2(v_hat));
            
            aux = myfilter.*fft2(u.*v);
            
            result = real(ifft2(aux));            
        case 'fourier2real'
            u_hat = u;
            v_hat = v;
            
            u = real(ifft2(u_hat));
            v = real(ifft2(v_hat));
            result = (u.*v);
            result(abs(result) < 1e-16) = 0;
    end

return