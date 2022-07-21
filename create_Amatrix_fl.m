function [A, x, M] = create_Amatrix_fl(N, L, s)

% N is the number of interior nodes (we consider zero Dirichlet boundary conditions)
% L is the length of the domain
% s is the fractional power

x = linspace(0,L,N+2);
h = x(2) - x(1);

A = zeros(N,N);
if s == 0.5
    for ii = 1:N
        for jj = 1:N
            if jj == ii
                A(ii,jj) = 8*log(2);
            elseif abs(jj-ii) == 1
                A(ii,jj) = 9*log(3) - 16*log(2);
            elseif abs(jj-ii) == 2
                A(ii,jj) = 56*log(2) - 36*log(3);
            else
                A(ii,jj) = -4*(jj-ii+1)^2*log(jj-ii+1) - 4*(jj-ii-1)^2*log(jj-ii-1) ...
                    +6*(jj-ii)^2*log(jj-ii) + (jj-ii+2)^2*log(jj-ii+2) + (jj-ii-2)^2*log(jj-ii-2);
            end
        end
    end
else
    for ii = 1:N
        for jj = 1:N
            if jj == ii
                A(ii,jj) = 2*(2^(3-2*s)-4);
            elseif abs(jj-ii) == 1
                A(ii,jj) = 3^(3-2*s)-2^(5-2*s)+7;
            else
                kk = abs(jj-ii);
                A(ii,jj) = 4*(kk+1)^(3-2*s) + 4*(kk-1)^(3-2*s) - 6*kk^(3-2*s) ...
                    - (kk+2)^(3-2*s) - (kk-2)^(3-2*s);
            end
        end
    end
    A = -h^(1-2*s) * A / (2*s*(1-2*s)*(1-s)*(3-2*s));
end
c1s = s*2^(2*s)*gamma((1+2*s)/2) / sqrt(pi) / gamma(1-s);
A = c1s*A;

% also compute the matrix for dynamic problems
M  = sparse(N,N);
for e = 1:N-1
    ind = e:e+1;
    M(ind,ind) = M(ind,ind) + h*1/6*[2 1; 1 2];
end