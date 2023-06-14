function [Im_new, Cs_new, C_new, Is_new, I_new, Rs_new, R_new] = protrude_extrapolate_diff_model3(Im, Cs, C, Is, I, Rs, R, F, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A)

    geom_p = geometry_factor(Im, g, k);
    A_cor = extapolate_C(F, Im);
    V = sum(Im(:));
    s = size(Im);
    vol_p = 0.5 + beta_V - (beta_V+gamma_V)*(1-1./(1+exp(alpha_V*(V-V0))));
    act_p = 0.5 - beta_A + (beta_A+gamma_A)*(1-1./(1+exp(alpha_A*(A_cor-A_act))));
    protr_p = vol_p.*act_p.*geom_p;
    ex = rand(s)<protr_p;
    Im_new = Im + ex;
    
    Rs_cor = extapolate_C(Rs, Im);
    Rs_ex = ex.*Rs_cor;
    Rs_new = Rs + Rs_ex;
    Rs_new = Rs_new - Im_new*sum(Rs_ex(:))/sum(Im_new(:));
    
    R_cor = extapolate_C(R, Im);
    R_ex = ex.*R_cor;
    R_new = R + R_ex;
    R_new = R_new - Im_new*sum(R_ex(:))/sum(Im_new(:));
    
    Rs_neg = Rs_new;
    Rs_neg(Rs_neg > 0) = 0;
    Rs_new(Rs_new < 0) = 0;
    
    R_neg = R_new;
    R_neg(R_neg > 0) = 0;
    R_new(R_new < 0) = 0;
    
    Rs_new = Rs_new + R_neg;
    R_new = R_new + Rs_neg;
    
    Is_cor = extapolate_C(Is, Im);
    Is_ex = ex.*Is_cor;
    Is_new = Is + Is_ex;
    Is_new = Is_new - Im_new*sum(Is_ex(:))/sum(Im_new(:));
    
    I_cor = extapolate_C(I, Im);
    I_ex = ex.*I_cor;
    I_new = I + I_ex;
    I_new = I_new - Im_new*sum(I_ex(:))/sum(Im_new(:));
    
    Is_neg = Is_new;
    Is_neg(Is_neg > 0) = 0;
    Is_new(Is_new < 0) = 0;
    
    I_neg = I_new;
    I_neg(I_neg > 0) = 0;
    I_new(I_new < 0) = 0;
    
    Is_new = Is_new + I_neg;
    I_new = I_new + Is_neg;
    
    Cs_cor = extapolate_C(Cs, Im);
    Cs_ex = ex.*Cs_cor;
    Cs_new = Cs + Cs_ex;
    Cs_new = Cs_new - Im_new*sum(Cs_ex(:))/sum(Im_new(:));
    
    C_cor = extapolate_C(C, Im);
    C_ex = ex.*C_cor;
    C_new = C + C_ex;
    C_new = C_new - Im_new*sum(C_ex(:))/sum(Im_new(:));
    
    Cs_neg = Cs_new;
    Cs_neg(Cs_neg > 0) = 0;
    Cs_new(Cs_new < 0) = 0;
    
    C_neg = C_new;
    C_neg(C_neg > 0) = 0;
    C_new(C_new < 0) = 0;
    
    Cs_new = Cs_new + C_neg;
    C_new = C_new + Cs_neg;
end
