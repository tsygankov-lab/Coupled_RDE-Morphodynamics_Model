function l = laplacian_DT(con, h, U)
    
    l = 0*U; 
    l(2:(end-1),2:(end-1)) = - U(2:(end-1),2:(end-1)).*con(2:(end-1),2:(end-1)) + ...
          con(3:end,2:(end-1)) + con(1:(end-2),2:(end-1)) + ...
          con(2:(end-1),3:end) + con(2:(end-1),1:(end-2));
    l = l/h^2;  
    
end