function AM = fl1D_parts(A, M)

P = M*(M+1)/2;
N = size(A,1);

L = N/M;

if ~(floor(L) == ceil(L))
    error('Wrong number of partitions')
end

AM{P} = [];
p=0;

D = diag(A) + sum(abs(A - diag(diag(A)))).';

% diagonal blocks
for ii = 1:M
    indii = (ii-1)*L + (1:L);
    p = p+1;
    AM{p} = sparse(N,N);
    AM{p}(indii, indii) = A(indii, indii);
    AM{p} = AM{p} - diag(diag(AM{p}));
    AM{p} = AM{p} - diag(sum(abs(AM{p})));
    
    Dii = zeros(1,N);
    Dii(indii) = D(indii);
    AM{p} = AM{p} + diag(Dii);
end

% off diagonal blocks
for ii = 1:M
    indii = (ii-1)*L + (1:L);
    for jj = 1:ii-1
        indjj = (jj-1)*L + (1:L);
        p = p+1;
        AM{p} = sparse(N,N);
        AM{p}(indii, indjj) = A(indii, indjj);
        AM{p}(indjj, indii) = A(indjj, indii);
        AM{p} = AM{p} - diag(sum(abs(AM{p})));
    end
end

% % check
% sumAM = sparse(N,N);
% for p = 1:P
%     sumAM = sumAM + AM{p};
% end
% norm(A - sumAM) / norm(A)