function out_im = subPixelateImg(in_im,p)
% OUT_IM = SUBPIXELATEIMG(IN_IM,P)
% Subpixelate an NxN image IN_IM to an (N.P)x(N.P) image OUT_IM.
% CA Smith 2019

    in_imsiz = size(in_im);
    ndim = length(in_imsiz);

    out_imsiz = p*in_imsiz;
    out_im = nan(out_imsiz);

    switch ndim

        case 3
            for i = 1:in_imsiz(1)
                for j = 1:in_imsiz(2)
                    for k = 1:in_imsiz(3)
                        out_im(p*(i-1)+1:p*(i-1)+p, ...
                            p*(j-1)+1:p*(j-1)+p, ...
                            p*(k-1)+1:p*(k-1)+p) = in_im(i,j,k);
                    end
                end
            end
        case 2
            for i = 1:in_imsiz(1)
                for j = 1:in_imsiz(2)
                    out_im(p*(i-1)+1:p*(i-1)+p, ...
                        p*(j-1)+1:p*(j-1)+p) = in_im(i,j);
                end
            end

    end

end