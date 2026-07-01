filename = 'optIC_IC_randfour_N1_128_N2_128_dt_1.0e-04_K_1.0e+03_ell1_1.02_ell2_1.02_T_1.00e-01_opt_1_tol_1e-06_cont_0_optT_1.00e-01.dat';
fid = fopen(filename,'rb');
A = fread(fid, Inf, 'double');
fclose(fid);
A = A(1:2:end) + 1i*A(2:2:end);
u_hat = reshape(A, 128, 128);
u = ifft2(u_hat);
imagesc(real(u)); colorbar
colormap("turbo")
v = fftshift(abs(fft2(u)).^2);
energyL2 = sum( v(:) )*(1.02*1.02*(2*pi)^2)/(128^2)^2;