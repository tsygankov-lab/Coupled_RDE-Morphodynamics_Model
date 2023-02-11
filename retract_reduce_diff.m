function [Im_new, Xa_new, Xi_new, Ba_new, Bi_new] = retract_reduce_diff(Im, Xa, Xi, Ba, Bi, F1, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A)
         
    geom_p = geometry_factor(~Im, g, k);
    A_cor = extapolate_C(F1, Im);
    V = sum(Im(:));
    s = size(Im);
    vol_p = 0.5 - beta_V + (beta_V+gamma_V)*(1-1./(1+exp(alpha_V*(V-V0))));
    act_p = 0.5 + beta_A - (beta_A+gamma_A)*(1-1./(1+exp(alpha_A*(A_cor-A_act))));
    retr_p = vol_p.*act_p.*geom_p;
    shr = rand(s)<retr_p;
    Im_new = Im - shr;
    
    Xa_shr = shr.*Xa;
    Xa_new = (Xa + sum(Xa_shr(:))/sum(Im_new(:))).*Im_new;
    
    Xi_shr = shr.*Xi;
    Xi_new = (Xi + sum(Xi_shr(:))/sum(Im_new(:))).*Im_new;
    
    Ba_shr = shr.*Ba;
    Ba_new = (Ba + sum(Ba_shr(:))/sum(Im_new(:))).*Im_new;
    
    Bi_shr = shr.*Bi;
    Bi_new = (Bi + sum(Bi_shr(:))/sum(Im_new(:))).*Im_new;
end
