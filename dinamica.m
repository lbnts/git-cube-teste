
function wp = dinamica (wibb, N, J)

wp = J \ (-Skew(wibb) * J * wibb + N);
