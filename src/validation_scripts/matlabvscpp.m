% Initial state
vA = v_step;   % reference state
vB = v_step;   % fused state

Nonlin_v0 = 0; % stage-1 uses beta_E(1)=0 anyway

% --- Reference MATLAB loop (consistent)
v_fwd2x = reshape( D1vec .* v_fwdstep, [ N , N_x2 ] );                       % f_x
v_fwd2y = reshape( D2vec .* v_fwdstep, [ N , N_x2 ] );                       % f_y
v_fwd2lap = reshape( Lap .* v_fwdstep, [ N , N_x2 ] );                       % lap(f)
% nonlinear terms and solution substeps
for k = 1:4
    v_step2 = reshape( vA, [ N , N_x2 ] );                               % z (adjoint variable)
    v_step2x = reshape( D1vec .* vA, [ N , N_x2 ] );                     % z_x
    v_step2y = reshape( D2vec .* vA, [ N , N_x2 ] );                     % z_y
    w1_r = multiply2D( v_step2x , v_fwd2x , 'fourier2real' );                   % z_x * f_x in physical space (pseudospectral)
    w1s_r = multiply2D( v_step2y , v_fwd2y , 'fourier2real' );                  % z_y * f_y in physical space (pseudospectral)
    w2_r = multiply2D( v_fwd2lap , v_step2 , 'fourier2real' );                  % lap(f) * z in physical space (pseudospectral)
    Nonlin_v1_r = (-1)*( w1_r + w1s_r ) - ( w2_r ) ;                            % - (z_x * f_x + z_y * f_y) - (lap(f) * z) in physical space
    Nonlin_v1 = multiply2D( fft2(Nonlin_v1_r) , fft2(Nonlin_v1_r) , 'dealias'); % dealias
    Nonlin_v1 = Nonlin_v1(:);

    vA = ( 1 + (dt * alpha_I(k) * Lin) ).^(-1) .* ...
        ( ( 1 - (dt * beta_I(k) * Lin) ) .* vA - ...
        (dt * alpha_E(k) * Nonlin_v1) - (dt * beta_E(k) * Nonlin_v0) );

    Nonlin_v0 = Nonlin_v1;
end

% --- Fused MEX
vB = ks_nonlinear_bwd( vB, v_fwdstep, D1vec, D2vec, Lap, N, Lin, dt, alpha_I, beta_I, alpha_E, beta_E);

% --- Compare (reference vs fused)
rel_err = norm(vA - vB) / norm(vA)