function x_tilde = WLLS(Anc,Theta_hat,Phi_hat)
%Plain LLS

[~,L] = size(Anc);

E = [1 0 0; 0 1 0];

A_mtx = [];

W = [];

b = [];

%Initialize w and W
w = ones(L,1);

W = kron(eye(2),diag(w));

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

x_bar = (pinv(A_mtx'*W'*W*A_mtx))*(A_mtx'*W'*b);

sum_dist = 0;

for i = 1:L
    sum_dist = sum_dist + norm(x_bar-Anc(1:3,i));
end

for i = 1:L
    w(i) = 1 - ((norm(x_bar-Anc(1:3,i)))/sum_dist);
end

W = kron(eye(2),diag(w));

x_tilde = (inv(A_mtx'*W*A_mtx))*(A_mtx'*W*b);



end

