function[zx]=zci(v)
% zci: zero-crossing indices

zx_temp = find(v(:).*circshift(v(:), [-1 0]) < 0);
    
    if isempty(zx_temp)
        zx=0;
        elseif  sign(v(1))*sign(v(end))<0
            zx=zx_temp(1:end-1); % rule out last element
        else
        zx=zx_temp;
    end
end