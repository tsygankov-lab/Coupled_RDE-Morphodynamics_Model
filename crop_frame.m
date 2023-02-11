function [Im, i_start, i_end, j_start, j_end] = crop_frame(Im_L, crop_d)

    i_start = find(sum(Im_L, 2),1);
    i_end = find(sum(Im_L, 2),1, 'last');

    j_start = find(sum(Im_L, 1),1);
    j_end = find(sum(Im_L, 1),1, 'last');

    i_start = i_start - crop_d;
    i_end = i_end + crop_d;
    j_start = j_start - crop_d;
    j_end = j_end + crop_d;

    Im = Im_L(i_start:i_end, j_start:j_end);

end