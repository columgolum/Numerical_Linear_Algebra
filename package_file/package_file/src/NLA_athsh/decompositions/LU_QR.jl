function LU(M :: AbstractMatrix{<:Number})
    U = Float64.(M)
    sz = size(U)
    tol = _zero_tol(U)  # anything this small is treated as an exact zero
    L = zeros(Float64, sz[1], sz[1])
    P = zeros(Float64, sz[1], sz[1])
  @inbounds @simd  for modi in 1:1:sz[1]
        L[modi, modi] = 1
        P[modi, modi] = 1
    end
    i = 1
    rank = 0
    h = 1 
    
    while i <= sz[1] && h <= sz[2] # traversing both rows and columns
        pivot_row = i
        pivot = abs(U[i, h])
        
        # Scan ALL rows below row i in the current column h
        for meow in i+1:1:sz[1]
            temp = abs(U[meow, h]) 
            if temp > pivot # partial pivoting is occuring here
                pivot = temp
                pivot_row = meow
            end
        end
        
        # If the largest element in this column is (numerically) 0 then skip to the next column
        # Because a number like 0 or near zero leads to large numbers ( when divided by 0 or near zero)
        # which leads to overflow and roundoff errors
        if pivot <= tol
            U[i:sz[1], h] .= 0.0  # clean the noise, so later "!= 0" checks are reliable
            h += 1
            continue
        end
        
        # pivot found hence increase the rank by 1
        rank += 1
        h_carry = h
        
        # rows of U swapped
        temp = U[i, :]
        U[i, :] = U[pivot_row, :]
        U[pivot_row, :] = temp
        
        # rows of P swapped
        temp = P[i, :]
        P[i, :] = P[pivot_row, :]
        P[pivot_row, :] = temp
        
        # rows of L swapped
        if i > 1
            temp = L[i, 1:i-1]
            L[i, 1:i-1] = L[pivot_row, 1:i-1]
            L[pivot_row, 1:i-1] = temp
        end 
        
        pivot_val = U[i, h_carry]
        @inbounds @simd for j in i+1:1:sz[1]
            mbp = U[j, h_carry] / pivot_val 
            L[j, i] = mbp
            @inbounds @simd for k = 1:1:sz[2]
                U[j, k] = U[j, k] - mbp * U[i, k]
            end
            U[j, h_carry] = 0.0  # mbp*pivot is not always EXACTLY the entry, so force the exact zero
        end
        
        i += 1
        h += 1
    end
    return P, L, U, rank
end

function Householder_matrix(v::Vector{<:Number})
    # Returns the n x n Householder matrix H = I - 2 w w' / (w' w) that maps v to -sign(v[1])*||v|| * e1.
    # The input v is not modified.
    n = length(v)
    H = zeros(Float64, n, n)
    for i in 1:n
        H[i, i] = 1.0
    end

    nrm = sqrt(sum(abs2, v))
    if nrm == 0                       # zero vector: nothing to reflect
        return H
    end

    w = Float64.(v)
    w[1] += v[1] >= 0 ? nrm : -nrm    # choose the sign that avoids cancellation
    wTw = sum(abs2, w)
    H .-= (2 / wTw) .* (w * transpose(w))
    return H
end

function QR(M_raw::AbstractMatrix{<:Number})
    # Returns Q (m x m, orthogonal) and R (m x n, upper triangular) with M_raw = Q * R.
    R = Float64.(M_raw)
    m, n = size(R)
    Q = zeros(Float64, m, m)
    for i in 1:m
        Q[i, i] = 1.0
    end

    # The last reflection on a single-entry column does nothing, so we stop at m-1.
    n_reflections = min(m - 1, n)
    for i in 1:n_reflections
        x = R[i:m, i]                 # copy of the sub-column that we want to zero out
        H = Householder_matrix(x)

        R[i:m, i:n] = H * R[i:m, i:n] # apply H to the active block of R
        Q[:, i:m] = Q[:, i:m] * H     # accumulate Q = H1 * H2 * ... * Hk
        R[i+1:m, i] .= 0.0            # clean up the entries below the diagonal
    end
    return Q, R
end

function solutions(U::Matrix{<:Number},b_raw::AbstractVecOrMat,r::Int64)
    rank = r
    sz = size(U)    

    particular_vector = zeros(Float64,sz[2],1)
    row_status = zeros(Int64,sz[1],2) 
    # rows of row_status = number of rows
    # first column tells wheather a row is pivot row(1) or free row(0) 
    # second column tells in which column that pivot is 
    column_status = zeros(Int64,sz[2],1)
    # columns of column_status = number of columns in A which is A|b = U
    # first row tells wheather a column is pivot column(1) or free column(0)
    
 @inbounds  for i in 1:1:sz[1]
    @inbounds     for j in 1:1:sz[2]
            if U[i,j] !=0
                row_status[i,1] = 1
                row_status[i,2] = j
                column_status[j] =1
                break
            
            end
        end
    end # this code block is used to find the pivot columns and pivot rows, and store them in row_status and column_status respectively

@inbounds for i in sz[1]:-1:1
        if row_status[i,1] ==0
             nothing
        else
            pivot_column = row_status[i,2]
            b_comp = b_raw[i]
          @inbounds @simd  for doraemon in pivot_column+1:1:sz[2]
                if column_status[doraemon] == 0
                    nothing
                else
                    b_comp -= U[i,doraemon]*particular_vector[doraemon]
                end
            end
            particular_vector[pivot_column]= b_comp/U[i,pivot_column]
        end
    end  # this code block is used to find the particular solution of the system of linear equations, and store it in particular_vector

    if rank == sz[2]
        return particular_vector
    else 

       linear_combination_of_vectors = Matrix{Float64}[]   
@inbounds for j in 1:1:sz[2]
            if column_status[j] == 0 #= We only care about free columns where column_status[j] == 0
                                         Because the pivot variabels of basis vector are (-1)*M[:,j] where j is a free column
                                        =#
            
            basis_vector = zeros(Float64,sz[2],1)
            basis_vector[j] = 1 # the variable of free column(free variable is set to 1 => this is linear algebra theory)
          @inbounds for i in sz[1]:-1:1
                if row_status[i,1] ==0
                    nothing
                else
                    pivot_column = row_status[i,2]
                    b_comp = 0.0
                    @inbounds @simd  for doraemon in pivot_column+1:1:sz[2]
                        b_comp -= U[i,doraemon]*basis_vector[doraemon]
                                     end
                    basis_vector[pivot_column]= b_comp/U[i,pivot_column]
                end
            end
            push!(linear_combination_of_vectors,basis_vector)
        end
    end # this code block is used to find the linear combination of vectors of the system of linear equations, and store it in linear_combination_of_vectors


return particular_vector, linear_combination_of_vectors   
        
end
end
