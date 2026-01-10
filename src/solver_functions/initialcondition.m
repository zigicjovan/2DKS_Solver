function u_0 = initialcondition(IC,utility1,x1,x2,N,N_x2,L_x1,L_x2,Klap)

    u_0 = 0;
    kabs = sqrt(Klap);              % absolute radii
    kmax = max(kabs(:));            % maximum radius
    switch IC 
        case 'optimized'
            u_0 = reshape( utility1, [ N , N_x2 ] );
        case 'noise'
            u_0 = 2*rand(N)-ones(N);
        case 'randfour'
            delta = 6/kmax;                 % AS width            
            A = exp(-delta * kabs);         % analytic envelope
            A(1,1) = 0;                     % zero mean (k=0)
            phi = 2*pi*rand( N , N_x2 );    % uniform random phase
            u_0 = A .* exp(1i*phi);         % Phase-randomized Fourier field
            %eta = (randn(size(A)) + 1i*randn(size(A))) / sqrt(2);
            %u_0 = A .* eta;         % Phase-randomized Fourier field
            u_0 = real(ifft2(u_0));         % enforce real field symmetry: uhat(-k) = conj(uhat(k))
        case 's'
            u_0 = sin( (x1 + x2) ) + sin( x1 ) + sin( x2 );
        case 's1'
            u_0 = sin( (L_x1*x1 + L_x2*x2) ) + sin( L_x1*x1 ) + sin( L_x2*x2 );
        case 'sc1'
            u_0 = cos( L_x1*x1 ) .* sin( L_x2*x2 );
        case 'shf'
            u_0 = sin( 1*(L_x1*x1 + L_x2*x2) ) + sin( kmax*L_x1*x1 ) + sin( kmax*L_x2*x2 );
        case 'shfn'
            u_0 = sin( 1*(L_x1*x1 + L_x2*x2) ) + sin( 0.1*kmax*L_x1*x1 ) + sin( 0.1*kmax*L_x2*x2 );
        case 'tg1'
            u_0 = sin( L_x1*x1 ) .* sin( L_x2*x2 );
        case 'tghf'
            u_0 = sin( L_x1*x1 ) .* sin( kmax*L_x2*x2 );
        case 'stg1'
            u_0 = sin( (L_x1*x1 + L_x2*x2) ) + sin( L_x1*x1 ) .* sin( L_x2*x2 );
        case 'stghf'
            u_0 = sin( (L_x1*x1 + L_x2*x2) ) +  sin( L_x1*x1 ) .* sin( kmax*L_x2*x2 );
        case 'mmhf'
            u_0 = sin( (L_x1*x1 + L_x2*x2) ) + sin( ceil(kmax/10)*L_x1*x1 ) + sin( ceil(kmax/5)*L_x2*x2 ) + ...
                     sin( ceil(1*kmax/2)*(L_x1*x1 + L_x2*x2) ) + sin( ceil(kmax/3)*L_x1*x1 ) + sin( ceil(1*kmax/4)*L_x2*x2 ) + ...
                     sin( ceil(2*kmax/3)*(L_x1*x1 + L_x2*x2) ) + sin( kmax*L_x1*x1 ) + sin( ceil(3*kmax/4)*L_x2*x2 ) ;
        case 'mn'
            u_0 = sin( 1.0*(L_x1*x1 + L_x2*x2) ) + sin( 5.0*L_x1*x1 ) + sin( 10.0*L_x2*x2 ) + ...
                     sin( 60.0*(L_x1*x1 + L_x2*x2) ) + sin( 50.0*L_x1*x1 ) + sin( 30.0*L_x2*x2 ) + ...
                     sin( 80.0*(L_x1*x1 + L_x2*x2) ) + sin( 150.0*L_x1*x1 ) + sin( 90.0*L_x2*x2 ) ;
    end
end