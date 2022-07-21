function [A, duration, connections] = create_AMrand(nx, ny, nz, dx, dy, dz, node_pos,M)

tic
nn = nx*ny*nz;
indi = zeros(M, 12*nn); indj = zeros(M, 12*nn); Aval = zeros(M, 12*nn);
connections = zeros(3*nn, 7); count = 0;
ll = zeros(M,1);
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz
            node0 = node_pos(ii,jj,kk);
            if ii+1<=nx
                node1 = node_pos(ii+1,jj,kk);
                m = randi([1,M]);
                ll(m) = ll(m)+1; ind = 4*(ll(m)-1)+(1:4);
                indi(m,ind) = [node0, node0, node1, node1];
                indj(m,ind) = [node0, node1, node0, node1];
                Aval(m,ind) = [-1, 1, 1, -1]/ dx^2;
                count = count + 1; connections(count, :) = [ii,jj,kk,ii+1,jj,kk, m];
            end
            if jj+1<=ny
                node2 = node_pos(ii,jj+1,kk);
                m = randi([1,M]);
                ll(m) = ll(m)+1; ind = 4*(ll(m)-1)+(1:4);
                indi(m,ind) = [node0, node0, node2, node2];
                indj(m,ind) = [node0, node2, node0, node2];
                Aval(m,ind) = [-1, 1, 1, -1]/ dy^2;
                count = count + 1; connections(count, :) = [ii,jj,kk,ii,jj+1,kk, m];
            end
            if kk+1<=nz
                node3 = node_pos(ii,jj,kk+1);
                m = randi([1,M]);
                ll(m) = ll(m)+1; ind = 4*(ll(m)-1)+(1:4);
                indi(m,ind) = [node0, node0, node3, node3];
                indj(m,ind) = [node0, node3, node0, node3];
                Aval(m,ind) = [-1, 1, 1, -1]/ dz^2;
                count = count + 1; connections(count, :) = [ii,jj,kk,ii,jj,kk+1, m];
            end
        end
    end
end

diagD = zeros(nn,1);
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz
            d = 1;
            if ii == 1 || ii == nx
                d = d*0.5;
            end
            if jj == 1 || jj == ny
                d = d*0.5;
            end
            if kk == 1 || kk == nz
                d = d*0.5;
            end
            node0 = node_pos(ii,jj,kk);
            diagD(node0) = d;
        end
    end
end

A = cell(1,M);
for m = 1:M
    indim = indi(m,1:4*ll(m)); indjm = indj(m,1:4*ll(m)); Avalm = Aval(m, 1:4*ll(m));
    A{m} = spdiags(diagD,0,nn,nn)\sparse(indim,indjm,Avalm,nn,nn);
end
connections = connections(1:count, :);

duration = toc;