function geom_p = geometry_factor(Im, g, k)

    % expand_weights v2
    %
    % changes: function has additional g (geometric power) parameter as input
    % this parameter means power of 1/distance, when weights are calculated to accaount for 
    % influence of neighboring pixels
    % the most logical value of this parameter is 2 (from physics, ~1/r^2)
    % additionaly, only r=1 is left as the raduis of interaction
    % 
    % function returns weigths for pixels where cell can expand
    %
    % cell is represented as 1s and 0s. transitions 1->0 or 0->1 are equal to
    % cell motion
    %
    % 1. transition is possible only in presence of direct neighbors (vertical
    % or horisontal), all diagonal transitions are prohibited (to avoid touching)
    %
    % 2. transition is prohibited if vertical direct neighbors are 1s and
    % horisontal direct neighbors are 0s and vice versa (to avoi bubbles formation)
    
    % R - right
    % L - left
    % U - up
    % D - down
    
    RIm = 0*Im;
    RIm(1:end, 2:end) = Im(1:end, 1:end-1);
    RIm = (RIm-Im)==1;
    
    LIm = 0*Im;
    LIm(1:end, 1:end-1) = Im(1:end, 2:end);
    LIm = (LIm-Im)==1;
    
    UIm = 0*Im;
    UIm(1:end-1, 1:end) = Im(2:end, 1:end);
    UIm = (UIm-Im)==1;
    
    DIm = 0*Im;
    DIm(2:end, 1:end) = Im(1:end-1, 1:end);
    DIm = (DIm-Im)==1;
    
    RUIm = 0*Im;
    RUIm(1:end-1, 2:end) = Im(2:end, 1:end-1);
    RUIm = (RUIm-Im)==1;
    
    RDIm = 0*Im;
    RDIm(2:end, 2:end) = Im(1:end-1, 1:end-1);
    RDIm = (RDIm-Im)==1;
    
    LUIm = 0*Im;
    LUIm(1:end-1, 1:end-1) = Im(2:end, 2:end);
    LUIm = (LUIm-Im)==1;
    
    LDIm = 0*Im;
    LDIm(2:end, 1:end-1) = Im(1:end-1, 2:end);
    LDIm = (LDIm-Im)==1;
    
    %diagonal elements
    RUIm_d = RUIm-RIm-UIm == 1;
    RDIm_d = RDIm-RIm-DIm == 1;
    LUIm_d = LUIm-LIm-UIm == 1;
    LDIm_d = LDIm-LIm-DIm == 1;
    diag = RUIm_d + RDIm_d + LUIm_d + LDIm_d;
    
    % avoid bubbles, prohibit transition if vertical neighbors are 1s and
    % horisontal are 0s and vice versa
    
    vert_hor = (((RIm + LIm) == 2) & ((UIm + DIm) == 0)) + (((RIm + LIm) == 0) & ((UIm + DIm) == 2));
    
    prohibit = (diag + vert_hor) >= 1;
   
    r = 1; %by default we consider that distance between pixels is 1 
    weights = ((RIm + LIm + UIm + DIm)/(r^g) + (RUIm + RDIm + LUIm + LDIm)/((r*sqrt(2))^g)).*(~prohibit)/(4/(r^g) + 4/((r*sqrt(2))^g));
    
    bw = weights>0;
    
    % eliminate only one element that can lead to diagonal touching
    % also eliminate direct touching in square 2x3 elements 
    bw(1,:) = 0;
    bw(end,:) = 0;
    bw(:,1) = 0;
    bw(:,end) = 0;
    ind = find(bw); % find non-zero elements
    rand_ind = randperm(length(ind)); % randomized order of elements
    for i=rand_ind
        [a,b] = ind2sub(size(bw), ind(i));
        
        %this helps to avoid simultaneous diagonal touching
        if ~(((~bw(a-1,b-1)) || Im(a, b-1) || Im(a-1, b)) && ...
            ((~bw(a+1, b-1)) || Im(a+1, b) || Im(a, b-1)) && ...
            ((~bw(a-1, b+1)) || Im(a-1, b) || Im(a, b+1)) && ...
            ((~bw(a+1, b+1)) || Im(a+1, b) || Im(a, b+1)))
        
            bw(a,b) = 0;
        end
        
        %this helps to avoid simultaneous direct touching
        if (bw(a+1, b) && ~Im(a+1, b+1) && ~Im(a+1, b-1) && ~Im(a, b+1) && ~Im(a, b-1)) || ...
           (bw(a-1, b) && ~Im(a-1, b+1) && ~Im(a-1, b-1) && ~Im(a, b+1) && ~Im(a, b-1)) || ...
           (bw(a, b+1) && ~Im(a+1, b+1) && ~Im(a-1, b+1) && ~Im(a+1, b) && ~Im(a-1, b)) || ...
           (bw(a, b-1) && ~Im(a+1, b-1) && ~Im(a-1, b-1) && ~Im(a+1, b) && ~Im(a-1, b))
            
            bw(a,b) = 0;
        end
    end
    
    geom_p = (weights.*bw).^k;
end
