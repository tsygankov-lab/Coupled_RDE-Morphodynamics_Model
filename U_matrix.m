function U = U_matrix(Im)

    U = 0*Im;
    U(2:(end-1),2:(end-1))=Im(3:end,2:(end-1))+Im(1:(end-2),2:(end-1))+...
        Im(2:(end-1),3:end)+Im(2:(end-1),1:(end-2));
    U(Im==0)=0;

end