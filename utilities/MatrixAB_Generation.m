function [A,B] = MatrixAB_Generation(graph_select,normalizeAB,n)

if graph_select == 1    % strongly connected
    if normalizeAB == 0
        A = -(n-1)*eye(n);
        B = ones(n,n)-eye(n);
    else
        A = (-(n-1)*eye(n))/(n-1);
        B = (ones(n,n)-eye(n))/(n-1);
    end
elseif graph_select == 2  % ring (bilateral connection)
    if normalizeAB == 0
        A = -(2)*eye(n);
        B = circshift(eye(n),1)+circshift(eye(n),-1);
    else
        A = (-(2)*eye(n))/2;
        B = (circshift(eye(n),1)+circshift(eye(n),-1))/2;
    end
elseif graph_select == 3  % ring (directed connection)
    if normalizeAB == 0
        A = -(1)*eye(n);
        B = circshift(eye(n),-1);
    else
        A = (-(1)*eye(n))/1;
        B = circshift(eye(n),-1)/1;
    end
else
    error('graph_select is incorrect')
end
end
