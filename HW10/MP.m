function [index_set, coef_set, residual] = MP(y, A, N)

    y = y(:);
    col_norm = sqrt(sum(abs(A).^2, 1)) + eps;
    A_n = A ./ col_norm;

    residual = y;
    index_set = zeros(1, N);

    for it = 1:N
        corr = A_n' * residual;

        if it > 1
            corr(index_set(1:it-1)) = 0;
        end

        [~, idx] = max(abs(corr));
        index_set(it) = idx;

        As = A(:, index_set(1:it));
        coef_ls = (As' * As) \ (As' * y);

        residual = y - As * coef_ls;
    end

    coef_set = coef_ls;
end


