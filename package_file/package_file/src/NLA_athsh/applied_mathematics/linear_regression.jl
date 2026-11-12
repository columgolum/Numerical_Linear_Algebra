function linear_regression(X_raw :: AbstractMatrix{<:Number}, y_raw :: AbstractVecOrMat{<:Number}; mode :: String = "QR", intercept :: Bool = true)
    #=
    Least squares : finds the coefficients b which make the sum of (y - X*b)^2 as small as possible.
    X = m x k matrix (every row is a sample, every column is a feature), y = m values.
    If intercept = true then a column of ones is added in front of X, so the first coefficient is the intercept.

    Returns b as a column matrix  (p x 1).   Three ways to find b :
    "normal" : solve the normal equations (X'X) b = X'y with our own solveLSE (LU)  -> fastest, but squares the condition number
    "QR"     : X = QR  =>  R b = Q'y (back substitution)                           -> accurate, the default
    "SVD"    : b = V * (1/Sigma) * U' * y  (pseudo-inverse)                        -> also works when the features are dependent
    =#
    X = Float64.(X_raw)
    y = Float64.(y_raw)[:]       # a column matrix (m x 1) becomes a plain vector
    m = size(X, 1)

    if length(y) != m
        throw(ArgumentError("X and y must have the same number of rows"))
    end

    if intercept
        X = hcat(ones(Float64, m), X)
    end
    p = size(X, 2)

    if mode == "normal"
        XT = transpose(X)
        result = solveLSE(hcat(XT * X, XT * y), "LU")
        if result[1] != "Unique solution"
            throw(ArgumentError("The features are linearly dependent, use mode = \"SVD\""))
        end
        return result[2]

    elseif mode == "QR"
        if m < p
            throw(ArgumentError("QR mode needs at least as many samples as coefficients, use mode = \"SVD\""))
        end
        Q, R = QR(X)
        c = transpose(Q) * y
        # R = [R1 ; 0] where R1 is p x p upper triangular, so  R1 b = (Q'y)[1:p]  is a (square) system
        result = solveLSE(hcat(R[1:p, 1:p], c[1:p]), "GE")
        if result[1] != "Unique solution"
            throw(ArgumentError("The features are linearly dependent, use mode = \"SVD\""))
        end
        return result[2]

    elseif mode == "SVD"
        U, Sigma, V, rank = SVD(X)
        cutoff = max(size(X)...) * eps(Float64) * Sigma[1, 1]
        b = zeros(Float64, p)
        for i in 1:1:min(size(X)...)
            if Sigma[i, i] > cutoff      # singular values which are (numerically) 0 are skipped
                b += (transpose(U[:, i]) * y / Sigma[i, i]) * V[:, i]
            end
        end
        return reshape(b, p, 1)

    else
        throw(ArgumentError("Supported modes are = normal, QR, SVD"))
    end
end

function predict_linear(X_raw :: AbstractMatrix{<:Number}, b :: AbstractVecOrMat{<:Number}; intercept :: Bool = true)
    # y_predicted = X * b  (with the column of ones in front of X if the model has an intercept)
    X = Float64.(X_raw)
    if intercept
        X = hcat(ones(Float64, size(X, 1)), X)
    end
    coefficients = Float64.(b)[:]
    if size(X, 2) != length(coefficients)
        throw(ArgumentError("Number of coefficients does not match the number of features"))
    end
    return reshape(X * coefficients, size(X, 1), 1)
end
