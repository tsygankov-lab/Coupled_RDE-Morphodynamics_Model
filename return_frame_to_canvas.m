function Im_L = return_frame_to_canvas(Im, S, i_start, i_end, j_start, j_end)

    Im_L = zeros(S);
    Im_L(i_start:i_end, j_start:j_end) = Im;

end