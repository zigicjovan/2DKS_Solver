function [error_2,error_inf,comptime] = temporalconvergence_2DKS(timestep, method, N, T, L_s1, L_s2)

    %%% (1) choose smallest step %%%
    dt = timestep(end); % size of smallest timestep

    %%% (2) load reference solution %%%

    [u_n, ~, time_n] = load_2DKSsolution('time_evolution', method, dt, T, N, L_s1, L_s2);
    u_Nsq = reshape( u_n(:,end) , N, N);
    error_2 = NaN( length(timestep) , 2 );
    error_inf = error_2;
    comptime = NaN( length(timestep) , 2 );

    %%% (3) error analysis %%%
    error_2( end , 1) = dt;
    error_inf( end , 1) = dt;
    error_2( end , 2) = 0;
    error_inf( end , 2) = 0;

    %%% (4) computational time %%%
    comptime( end , 1) = dt;
    comptime( end , 2) = time_n;

    j = 1; % iteration counter for saved measures

    for i = 1 : length(timestep)-1 % plot linf, l2 errors againt dt_small

        %%% (1) choose larger step %%%
        dt = timestep(i); % size of comparison timestep

        %%% (2) load solution %%%
        [u_n, ~, time_n] = load_2DKSsolution('time_evolution', method, dt, T, N, L_s1, L_s2);
        u_i = reshape( u_n(:,end) , N, N);

        %%% (3) error analysis %%%
        error_2( j , 1) = dt;
        error_inf( j , 1) = dt;
        error_2( j , 2) = norm( u_i - u_Nsq , 2) / norm(u_Nsq , 2);
        error_inf( j , 2) = norm( u_i - u_Nsq , inf) / norm(u_Nsq , inf);

        %%% (4) computational time %%%
        comptime( j , 1) = dt;
        comptime( j , 2) = time_n;

        j = j + 1;

    end

end