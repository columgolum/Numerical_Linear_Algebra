function Identity_matrix(n::Int)
    I_n = zeros(Float64, n, n)
    for i in 1:1:n
        I_n[i, i] = 1.0
    end
    return I_n
end

function Givens_rotation(a::Float64, b::Float64)
    # Returns c and s such that    [ c  s ] [a]   [r]
    #                              [-s  c ] [b] = [0]     where r = sqrt(a^2 + b^2)
    # A Givens rotation is a 2x2 rotation which kills ONE entry (here b) without touching the length of the vector
    r = hypot(a, b)
    if r == 0
        return 1.0, 0.0       # both entries are already 0, nothing to rotate
    end
    return a / r, b / r
end

function rotate_columns!(X::Matrix{Float64}, p::Int, q::Int, c::Float64, s::Float64)
    # X = X * R  where R is a rotation in the (p,q) plane :
    # new column p =  c * (column p) + s * (column q)
    # new column q = -s * (column p) + c * (column q)
    @inbounds for i in 1:1:size(X, 1)
        xp = X[i, p]
        xq = X[i, q]
        X[i, p] = c * xp + s * xq
        X[i, q] = -s * xp + c * xq
    end
end

function rotate_rows!(X::Matrix{Float64}, p::Int, q::Int, c::Float64, s::Float64)
    # X = L * X  where L is a rotation in the (p,q) plane :
    # new row p =  c * (row p) + s * (row q)
    # new row q = -s * (row p) + c * (row q)
    @inbounds for j in 1:1:size(X, 2)
        xp = X[p, j]
        xq = X[q, j]
        X[p, j] = c * xp + s * xq
        X[q, j] = -s * xp + c * xq
    end
end

function SVD(M_raw :: AbstractMatrix{<:Number}; tol=1e-14, max_iter=500)
    #=
    Returns U (m x m, orthogonal), Sigma (m x n, diagonal), V (n x n, orthogonal) and rank with
    M_raw = U * Sigma * transpose(V). The singular values on the diagonal of Sigma are sorted from the biggest to the smallest.

    Golub-Kahan algorithm:
    1) Bidiagonalisation : Householder matrices from the LEFT and the RIGHT turn A into an upper bidiagonal matrix B
                           (only the diagonal and the super-diagonal are non-zero)   A = U * B * V'
    2) Implicit shifted QR steps on B : a Givens rotation from the right creates a BULGE (an unwanted non-zero entry)
       and then Givens rotations from the left/right chase the bulge down the matrix until it falls out.
       Every step makes the super-diagonal smaller, until it becomes 0 and the diagonal is the singular values.
    =#
    A = Float64.(M_raw)
    m, n = size(A)

    # We work with m >= n, if the matrix is wide then transpose it:  A' = U S V'  =>  A = V S' U'
    if m < n
        U_t, Sigma_t, V_t, rank_t = SVD(Matrix(transpose(A)); tol = tol, max_iter = max_iter)
        return V_t, Matrix(transpose(Sigma_t)), U_t, rank_t
    end

    # Scale the matrix so that the biggest entry is 1, this keeps squares like d^2 away from overflow/underflow
    # (the singular values are multiplied by "big" again at the end)
    big = maximum(abs, A)
    if big == 0
        return Identity_matrix(m), zeros(Float64, m, n), Identity_matrix(n), 0
    end
    B = A ./ big
    U = Identity_matrix(m)
    V = Identity_matrix(n)

    # ------------------------------------------------------------------------------------------------------
    # STEP 1 : bidiagonalisation
    for i in 1:1:n
        # Left Householder : makes the entries below the diagonal in column i zero
        if i < m   # for a single entry column there is nothing to zero out
            x = B[i:m, i]
            H = Householder_matrix(x)
            B[i:m, i:n] = H * B[i:m, i:n]
            U[:, i:m] = U[:, i:m] * H
            B[i+1:m, i] .= 0.0
        end

        # Right Householder : makes the entries to the right of the super-diagonal in row i zero
        if i <= n - 2
            x = B[i, i+1:n]
            H = Householder_matrix(x)
            B[i:m, i+1:n] = B[i:m, i+1:n] * H
            V[:, i+1:n] = V[:, i+1:n] * H
            B[i, i+2:n] .= 0.0
        end
    end

    # ------------------------------------------------------------------------------------------------------
    # STEP 2 : Golub-Kahan SVD step (Givens rotations chase the bulge) on the n x n bidiagonal part
    tiny = eps(Float64) * maximum(abs, B)   # a diagonal entry smaller than this is treated as an exact 0
    hi = n
    iter = 0

    while hi > 1
        # deflation : if a super-diagonal entry is negligible compared to its neighbours on the diagonal then make it 0
        for i in 1:1:hi-1
            if abs(B[i, i+1]) <= tol * (abs(B[i, i]) + abs(B[i+1, i+1]))
                B[i, i+1] = 0.0
            end
        end

        # last super-diagonal entry is 0 => B[hi,hi] is a singular value, chop the last row and column
        if B[hi-1, hi] == 0.0
            hi -= 1
            iter = 0
            continue
        end

        # find the active block lo:hi, in which all the super-diagonal entries are non-zero
        lo = hi - 1
        while lo > 1 && B[lo-1, lo] != 0.0
            lo -= 1
        end

        # Special case 1 : a diagonal entry B[k,k] = 0 (k is not the last one).
        # Row k is then just one entry B[k,k+1], we rotate it away using the rows below, the bulge moves to the right
        # along row k, until it falls out of the matrix. After that row k is entirely 0 and the block splits in two.
        zero_found = false
        for k in lo:1:hi-1
            if abs(B[k, k]) <= tiny
                B[k, k] = 0.0
                for j in k+1:1:hi
                    c, s = Givens_rotation(B[j, j], B[k, j])
                    rotate_rows!(B, j, k, c, s)
                    rotate_columns!(U, j, k, c, s)
                    B[k, j] = 0.0
                end
                zero_found = true
                break
            end
        end
        if zero_found
            continue
        end

        # Special case 2 : the last diagonal entry B[hi,hi] = 0.
        # Same idea but with columns : the bulge moves UP along column hi
        if abs(B[hi, hi]) <= tiny
            B[hi, hi] = 0.0
            for j in hi-1:-1:lo
                c, s = Givens_rotation(B[j, j], B[j, hi])
                rotate_columns!(B, j, hi, c, s)
                rotate_columns!(V, j, hi, c, s)
                B[j, hi] = 0.0
            end
            continue
        end

        if iter >= max_iter
            error("SVD failed to converge at submatrix size $hi")
        end

        # Wilkinson shift : the eigenvalue of the bottom right 2x2 block of B'B which is closer to its last entry
        a = B[hi-1, hi-1]^2 + (hi - 1 > lo ? B[hi-2, hi-1]^2 : 0.0)
        b = B[hi-1, hi-1] * B[hi-1, hi]
        d = B[hi, hi]^2 + B[hi-1, hi]^2
        delta = (a - d) / 2
        sgn = delta >= 0 ? 1.0 : -1.0
        denom = abs(delta) + sqrt(delta^2 + b^2)
        shift = denom == 0 ? d : d - sgn * b^2 / denom

        # first rotation (from the right) is decided by the first column of (B'B - shift*I), it creates the bulge at B[lo+1,lo]
        y = B[lo, lo]^2 - shift
        z = B[lo, lo] * B[lo, lo+1]

        for k in lo:1:hi-1
            if k > lo
                # the bulge sits at B[k-1,k+1], kill it with a rotation from the right
                y = B[k-1, k]
                z = B[k-1, k+1]
            end
            c, s = Givens_rotation(y, z)
            rotate_columns!(B, k, k+1, c, s)
            rotate_columns!(V, k, k+1, c, s)
            if k > lo
                B[k-1, k+1] = 0.0
            end

            # that rotation created a bulge at B[k+1,k], kill it with a rotation from the left
            c, s = Givens_rotation(B[k, k], B[k+1, k])
            rotate_rows!(B, k, k+1, c, s)
            rotate_columns!(U, k, k+1, c, s)
            B[k+1, k] = 0.0
            # this one creates the next bulge at B[k,k+2] (if there is a k+2), which is handled in the next loop
        end
        iter += 1
    end

    # ------------------------------------------------------------------------------------------------------
    # singular values must be >= 0 : if one is negative then flip the sign of the matching column of V
    d = zeros(Float64, n)
    for i in 1:1:n
        d[i] = B[i, i]
        if d[i] < 0
            d[i] = -d[i]
            V[:, i] = -V[:, i]
        end
    end

    # sort from biggest to smallest, and sort the columns of U and V in the same way
    order = sortperm(d, rev = true)
    d = d[order]
    U[:, 1:n] = U[:, order]
    V = V[:, order]

    Sigma = zeros(Float64, m, n)
    for i in 1:1:n
        Sigma[i, i] = d[i] * big
    end

    # rank = number of singular values which are not (numerically) 0
    # (the standard cutoff for SVD : biggest singular value * max(m,n) * machine epsilon)
    cutoff = max(m, n) * eps(Float64) * Sigma[1, 1]
    rank = 0
    for i in 1:1:n
        if Sigma[i, i] > cutoff
            rank += 1
        end
    end
    return U, Sigma, V, rank
end
