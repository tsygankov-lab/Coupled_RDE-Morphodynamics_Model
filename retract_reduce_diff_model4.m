function [Im_new, Rs_new, R_new, Is_new, I_new, Cs_new, C_new] = retract_reduce_diff_model4(Im0, Im, Rs, R, Is, I, Cs, C, F, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A)
         
    geom_p = geometry_factor(~Im, g, k);
    A_cor = extapolate_C(F, Im);
    V = sum(Im(:));
    s = size(Im);
    vol_p = 0.5 - beta_V + (beta_V+gamma_V)*(1-1./(1+exp(alpha_V*(V-V0))));
    act_p = 0.5 + beta_A - (beta_A+gamma_A)*(1-1./(1+exp(alpha_A*(A_cor-A_act))));
    retr_p = vol_p.*act_p.*geom_p;
    shr = rand(s)<retr_p;
    shr = shr.*(~Im0);
    Im_new = Im - shr;
    
    Rs_shr = shr.*Rs;
    Rs_new = (Rs + sum(Rs_shr(:))/sum(Im_new(:))).*Im_new;
    
    R_shr = shr.*R;
    R_new = (R + sum(R_shr(:))/sum(Im_new(:))).*Im_new;
    
    Is_shr = shr.*Is;
    Is_new = (Is + sum(Is_shr(:))/sum(Im_new(:))).*Im_new;
    
    I_shr = shr.*I;
    I_new = (I + sum(I_shr(:))/sum(Im_new(:))).*Im_new;
    
    Cs_shr = shr.*Cs;
    Cs_new = (Cs + sum(Cs_shr(:))/sum(Im_new(:))).*Im_new;
    
    C_shr = shr.*C;
    C_new = (C + sum(C_shr(:))/sum(Im_new(:))).*Im_new;
end
