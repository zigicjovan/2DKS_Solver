function result = multiply2D(u,v,space)

    switch space         
        case 'dealias'

            Nx = size(u,1);
            Ny = size(u,2);
        
            % Wavenumbers in FFT ordering
            kx = [0:Nx/2 -Nx/2+1:-1];
            ky = [0:Ny/2 -Ny/2+1:-1];
            [KX, KY] = meshgrid(kx, ky);
        
            % 2/3-rule cutoff
            kx_cut = (2/3) * (Nx/2);
            ky_cut = (2/3) * (Ny/2);
        
            mask = (abs(KX) <= kx_cut) & (abs(KY) <= ky_cut);
        
            result = u .* mask;
            result(abs(result) < 1e-16) = 0;           
        case 'fourier2real'
            u_hat = u;
            v_hat = v;
            
            u = real(ifft2(u_hat));
            v = real(ifft2(v_hat));
            result = (u.*v);
            result(abs(result) < 1e-16) = 0;
    end

return