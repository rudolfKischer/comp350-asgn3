function [U, L] = elimination(U, k, L)
    n = size(U, 1);
    for i = k + 1:n
        [U, L] = eliminate(U, k, i, L);
    end
end


function [U, L] = eliminate(U, k, i, L)
    [n, ~] = size(U);
    m_ik = U(i, k) / U(k, k);
    if U(k, k) < 0
        m_ik = m_ik * -1;
    end
    L(i, k) = m_ik;
    U(i, k:n) = U(i, k:n) - m_ik * U(k, k:n);
end