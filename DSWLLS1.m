function x_tilde = DSWLLS1(Anc,Theta_hat,Phi_hat,N)
%Algorithm 1

[~,L] = size(Anc);

E = [1 0 0; 0 1 0];

x_tilde_N = WLLS(Anc,Theta_hat,Phi_hat);

x_tilde = x_tilde_N;

A_mtx = [];

b = [];

k = [0;0;1];

u_hat_mtx = zeros(3, L);
c_hat_mtx = zeros(3, L);

for i = 1:L
    u_hat_mtx(1:3,i) = [cos(Phi_hat(i))*cos(Theta_hat(i));cos(Phi_hat(i))*sin(Theta_hat(i));sin(Phi_hat(i))];
    c_hat_mtx(1:3,i) = [-sin(Theta_hat(i));cos(Theta_hat(i));0];
    A_mtx = [A_mtx;c_hat_mtx(1:3,i)'];
    b = [b;c_hat_mtx(1:3,i)'*Anc(1:3,i)];
end

for i = 1:L
    A_mtx = [A_mtx;(k-u_hat_mtx(1:3,i)*sin(Phi_hat(i)))'];
    b = [b;(k-u_hat_mtx(1:3,i)*sin(Phi_hat(i)))'*Anc(1:3,i)];
end

deltaN = (A_mtx*x_tilde_N - b)'*(A_mtx*x_tilde_N - b)/L;

vec_to_be_taken = 1:L;

C = nchoosek(vec_to_be_taken,N);

cnt = 0;

for i = 1:(factorial(L))/(factorial(N)*factorial(L-N))
    
    cnt = cnt + 1;
    
    Anc_new = [];
    
    Theta_hat_new = [];
    
    Phi_hat_new = [];
    
    A_mtx_new = [];
    
    b_new = [];
    
    u_hat_mtx_new = zeros(3, N);
    c_hat_mtx_new = zeros(3, N);
    
    for N_idx = 1:N
    
        Anc_new = [Anc_new, Anc(:,C(cnt,N_idx))];
        
        Theta_hat_new = [Theta_hat_new; Theta_hat(C(cnt,N_idx))];
        
        Phi_hat_new = [Phi_hat_new; Phi_hat(C(cnt,N_idx))];
    
    end
    
    x_est_new = WLLS(Anc_new,Theta_hat_new,Phi_hat_new);
    
    for N_idx = 1:N
        
        u_hat_mtx_new(1:3,N_idx) = [cos(Phi_hat_new(N_idx))*cos(Theta_hat_new(N_idx));cos(Phi_hat_new(N_idx))*sin(Theta_hat_new(N_idx));sin(Phi_hat_new(N_idx))];
        c_hat_mtx_new(1:3,N_idx) = [-sin(Theta_hat_new(N_idx));cos(Theta_hat_new(N_idx));0];
        A_mtx_new = [A_mtx_new;c_hat_mtx_new(1:3,N_idx)'];
        b_new = [b_new;c_hat_mtx_new(1:3,N_idx)'*Anc_new(1:3,N_idx)];
        
    end
    
    for N_idx = 1:N
        A_mtx_new = [A_mtx_new;(k-u_hat_mtx_new(1:3,N_idx)*sin(Phi_hat_new(N_idx)))'];
        b_new = [b_new;(k-u_hat_mtx_new(1:3,N_idx)*sin(Phi_hat_new(N_idx)))'*Anc_new(1:3,N_idx)];
    end
    
    if (A_mtx_new*x_est_new - b_new)'*(A_mtx_new*x_est_new - b_new)/N < deltaN 
        
        deltaN = (A_mtx_new*x_est_new - b_new)'*(A_mtx_new*x_est_new - b_new)/N;
        
        x_tilde = x_est_new;
        
    end
    
end




end

