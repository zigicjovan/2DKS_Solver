function checkCFL(dx,dt,u)
    C = 0.9;
    u_Linf = max(max(abs(u)));
    cflstep = C*dx/(max(1e-16,u_Linf));
    if dt > cflstep
        error(['Warning: CFL condition not upheld for solution coordinate = ' num2str(u_Linf) ' and dt = ' num2str(dt)]);
    end
end
