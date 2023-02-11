function outline = outline_8p(Im)

    outline = zeros(size(Im));
    outline(2:end-1,2:end-1) = Im(1:end-2,2:end-1) + Im(3:end,2:end-1) + ...
        Im(2:end-1,1:end-2) + Im(2:end-1,3:end) + ...
        Im(1:end-2,1:end-2) + Im(3:end,3:end) + ...
        Im(3:end,1:end-2) + Im(1:end-2,3:end);
    %outline = outline - 8*Im;
    outline(Im==1) = 0;
    outline = double(outline>0);

end