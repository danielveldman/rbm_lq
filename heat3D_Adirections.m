function [A1, A2, A3, duration] = heat3D_Adirections(nx, ny, nz, dx, dy, dz, node_pos)

tic
nn = nx*ny*nz;

%% less efficient but more insightful implementation
% A = sparse(nn,nn);
% for ii = 1:nx
%     for jj = 1:ny
%         for kk = 1:nz
%             node0 = node_pos(ii,jj,kk);
%             if ii+1<=nx
%                 node1 = node_pos(ii+1,jj,kk); ind = [node0, node1];
%                 A(ind,ind) = A(ind,ind) + [-1 1; 1 -1]/dx^2;
%             end
%             if jj+1<=ny
%                 node2 = node_pos(ii,jj+1,kk); ind = [node0, node2];
%                 A(ind,ind) = A(ind,ind) + [-1 1; 1 -1]/dy^2;
%             end
%             if kk+1<=nz
%                 node3 = node_pos(ii,jj,kk+1); ind = [node0, node3];
%                 A(ind,ind) = A(ind,ind) + [-1 1; 1 -1]/dz^2;
%             end
%         end
%     end
% end
% toc


indi = zeros(1, 4*nn); indj = zeros(1, 4*nn); Aval = zeros(1, 4*nn);
ll = 0;
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz
            node0 = node_pos(ii,jj,kk);
            if ii+1<=nx
                node1 = node_pos(ii+1,jj,kk);
                ll = ll+1; ind = 4*(ll-1)+(1:4);
                indi(ind) = [node0, node0, node1, node1];
                indj(ind) = [node0, node1, node0, node1];
                Aval(ind) = [-1, 1, 1, -1]/ dx^2;
            end
        end
    end
end
indi = indi(1:4*ll); indj = indj(1:4*ll); Aval = Aval(1:4*ll);
A1 = sparse(indi,indj,Aval);

indi = zeros(1, 4*nn); indj = zeros(1, 4*nn); Aval = zeros(1, 4*nn);
ll = 0;
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz
            node0 = node_pos(ii,jj,kk);
            if jj+1<=ny
                node2 = node_pos(ii,jj+1,kk);
                ll = ll+1; ind = 4*(ll-1)+(1:4);
                indi(ind) = [node0, node0, node2, node2];
                indj(ind) = [node0, node2, node0, node2];
                Aval(ind) = [-1, 1, 1, -1]/ dy^2;
            end
        end
    end
end
indi = indi(1:4*ll); indj = indj(1:4*ll); Aval = Aval(1:4*ll);
A2 = sparse(indi,indj,Aval);

indi = zeros(1, 4*nn); indj = zeros(1, 4*nn); Aval = zeros(1, 4*nn);
ll = 0;
for ii = 1:nx
    for jj = 1:ny
        for kk = 1:nz
            node0 = node_pos(ii,jj,kk);
            if kk+1<=nz
                node3 = node_pos(ii,jj,kk+1);
                ll = ll+1; ind = 4*(ll-1)+(1:4);
                indi(ind) = [node0, node0, node3, node3];
                indj(ind) = [node0, node3, node0, node3];
                Aval(ind) = [-1, 1, 1, -1]/ dz^2;
            end
        end
    end
end
indi = indi(1:4*ll); indj = indj(1:4*ll); Aval = Aval(1:4*ll);
A3 = sparse(indi,indj,Aval);

duration = toc;