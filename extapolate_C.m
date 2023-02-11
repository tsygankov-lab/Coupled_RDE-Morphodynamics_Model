function A_cor = extapolate_C(A,Im)

    A_cor = 0*A;
    A_cor(2:end-1,2:end-1) = (A(2:end-1,1:end-2) + A(2:end-1,3:end) + ...
        A(1:end-2,1:end-2) + A(1:end-2,2:end-1) + A(1:end-2,3:end) + ...
        A(3:end,1:end-2) + A(3:end,2:end-1) + A(3:end,3:end));
    
    %the folowing steps are made to avoide 1/Inf
    y = zeros(size(Im));
    y(2:end-1,2:end-1) = (Im(2:end-1,1:end-2) + Im(2:end-1,3:end) + ...
        Im(1:end-2,1:end-2) + Im(1:end-2,2:end-1) + Im(1:end-2,3:end) + ...
        Im(3:end,1:end-2) + Im(3:end,2:end-1) + Im(3:end,3:end));
    yy = 1./y;
    yy(yy==Inf) = 0;
    A_cor = A_cor.*yy;

end