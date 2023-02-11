function [Im_new, Xa_new, Xi_new, Ba_new, Bi_new] = protrude(Im0, Im, Xa, Xi, Ba, Bi, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A);
    geom_p = geometry_factor(Im, g, k);
    A_cor = extapolate_C(Xa, Im);
    V = sum(Im(:));
    s = size(Im);
    vol_p = 0.5 + beta_V - (beta_V+gamma_V)*(1-1./(1+exp(alpha_V*(V-V0))));
    act_p = 0.5 - beta_A + (beta_A+gamma_A)*(1-1./(1+exp(alpha_A*(A_cor-A_act))));
    protr_p = vol_p.*act_p.*geom_p;
    ex = rand(s)<protr_p;
    ex = ex.*(~Im0);
    Im_new = Im + ex;
    
    Xa_ex = ex.*A_cor;
    Xa_new = Xa + Xa_ex;
    
    Xi_cor = extapolate_C(Xi, Im);
    Xi_ex = ex.*Xi_cor;
    Xi_new = Xi + Xi_ex;
    Xi_new = Xi_new - Im_new*sum(Xa_ex(:)+Xi_ex(:))/sum(Im_new(:));
    
    Ba_cor = extapolate_C(Ba, Im);
    Ba_ex = ex.*Ba_cor;
    Ba_new = Ba + Ba_ex;
    
    Bi_cor = extapolate_C(Bi, Im);
    Bi_ex = ex.*Bi_cor;
    Bi_new = Bi + Bi_ex;
    Bi_new = Bi_new - Im_new*sum(Ba_ex(:)+Bi_ex(:))/sum(Im_new(:));
end
