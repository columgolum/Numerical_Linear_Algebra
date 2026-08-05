function Hessenberg_form(M_raw :: AbstractMatrix{<:Number})
    #= 
    Converting a matrix to upper Hessenberg form makes the QR shift algorithm  faster
    We do this by applying  Householder matrix  from both sides (Similarity transform to preserve eigenvalues)
    =#
    H = Float64.(M_raw)
    sz = size(H)
    
    if sz[1] > 2
     @inbounds @simd for i in 1:1:(sz[1]-2)
            # we only go up to sz[1]-2 because the last two columns don't need to be changed
            column_vec = @view H[i+1:sz[1], i] # selecting the sub-column for householder
            
            #------------------------------------------------------------------------------------------------------
            hehe = Householder_matrix(copy(column_vec)) 
            
            # Left multiply (applies to rows i+1 to end)
            H[i+1:sz[1], i:sz[2]] = hehe * H[i+1:sz[1], i:sz[2]] 
            
            # Right multiply (applies to columns i+1 to end)
            # Since hehe is symmetric, hehe = transpose(hehe). 
            # A similarity transform is Q * A * Q^T, so we multiply on the right as well, this preserves eigenvalues
            H[1:sz[1], i+1:sz[2]] = H[1:sz[1], i+1:sz[2]] * hehe
            
            #---------------------------------------------------------------------------------------------
            H[i+2:sz[1], i] .= 0.0 # making elements below SUB-diagonal zero
        end
    end
    return H
end

function eigenvalues(M_raw :: AbstractMatrix{<:Number}; tol=1e-12, max_iter=1000)
    A = Float64.(M_raw)
    sz = size(A)
    n = sz[1]
    
    # 1. Reduce matrix to upper Hessenberg form using Householder logic
    H = Hessenberg_form(A)
    low = n
    
    eigenvals = zeros(ComplexF64, n) 
    iter = 0
    
   while low > 1
        q = low - 1
        
        # checking deflation ( to avoid excessive computation)
        scale1 = abs(H[q, q]) + abs(H[low, low])
        if abs(H[low, q]) <= tol * max(scale1, 1.0) 
            H[low, q] = 0.0    # if a_i+1,i <= e*(aii+ai+1,i+1) then set ai+1,i = 0 
                                # and eigenvalue found = ai+1,i+1,
                                # chop the (i+1)th row and column 
            low -= 1        
            iter = 0
            continue
        end
    
        if low == 2
            # The entire remaining matrix is just a 2x2 block. 
            # We don't need to iterate further, the extraction loop will solve it exactly
            break 
        end
        
        p = low - 2
        scale2 = abs(H[p, p]) + abs(H[q, q])
        if abs(H[q, p]) <= tol * max(scale2, 1.0)
            H[q, p] = 0.0 # set  the sub-diagonal element above the 2x2 block to zero
            low -= 2      # shrink our effective matrix size by 2!
            iter = 0
            continue
        end
        
        if iter >= max_iter
            error(" QR algorithm failed to converge at submatrix size $low")
        end

         # Wilkinson shift calculation using the bottom-right 2x2 block
        a = H[q, q]
        b = H[q, low]
        c = H[low, q]
        d = H[low, low]
        
        trace_val = a + d
        det_val = a * d - b * c
        discrim = trace_val^2 - 4 * det_val
        
        # Keep the shift real,  If we used a complex shift here, we would have to 
        # convert the entire active_sub matrix to ComplexF64 for the QR step.
        if discrim >= 0
            sq = discrim^0.5
        else
            sq = 0.0 
        end
        
        lamda1 = (trace_val + sq) / 2.0
        lamda2 = (trace_val - sq) / 2.0
        
        # shift is whichever eigenvalue of the 2x2 block is closer to the bottom right element 'd'
        if abs(lamda1 - d) < abs(lamda2 - d)
            shift = lamda1
        else 
            shift =lamda2
        end
        
        # Construct shifted active submatrix => H_active - shift * I
        active_sub = copy(H[1:low, 1:low])
        
     @inbounds @simd for i in 1:1:low
            active_sub[i, i] -= shift
        end
    
        # decompose the active matrix 
        Q_sub, R_sub = QR(active_sub)
        
        # Compute R * Q + shift * I = new matrix
        next_sub = R_sub * Q_sub
        
     @inbounds @simd for i in 1:1:low
            next_sub[i, i] += shift
        end
        #------------------------------------------------------------------------------------------------------

         H[1:low, 1:low] = next_sub
        iter += 1
    end
    
    # Extract eigenvalues from the diagonal and the 2x2 blocks
    i = 1
    while i <= n
        if i < n && abs(H[i+1, i]) > tol
            a = H[i, i]
            b = H[i, i+1]
            c = H[i+1, i]
            d = H[i+1, i+1]
            
            trace_val = a + d
            det_val = a * d - b * c
            discrim = trace_val^2 - 4 * det_val
            
            if discrim >= 0
                sq = sqrt(discrim)
                eigenvals[i] = (trace_val + sq) / 2.0
                eigenvals[i+1] = (trace_val - sq) / 2.0
            else
                # Cast discrim to Complex so sqrt() knows to return an imaginary number
                sq_c = sqrt(Complex(discrim)) 
                eigenvals[i] = (trace_val + sq_c) / 2.0
                eigenvals[i+1] = (trace_val - sq_c) / 2.0
            end
            i += 2
        else
            eigenvals[i] = H[i, i]
            i += 1
        end
    end
    
    return eigenvals
end

function ref_eigen(M::Matrix{<:Number}; tol::Float64=1e-10)
    #=
    this code is very similar to ref in RDIS
    =#
    U = ComplexF64.(M) # Forces the matrix to allow complex numbers
    sz = size(U)
    i = 1
    rank = 0
    h = 1 
    
    while i <= sz[1] && h <= sz[2] 
        pivot_row = i
        pivot = abs(U[i, h])
        
        # Scan ALL rows below row i in the current column h
        segment = @view U[i+1:sz[1], h]
      @inbounds @simd for meow in 1:1:sz[1]-i
            temp = abs(segment[meow]) 
            if temp > pivot # partial pivoting is occurring here
                pivot = temp
                pivot_row = meow + i
            end
        end
        
        # Using  tolerance instead of strict == 0  as compared to REF in RDIS file
        if pivot <= tol
            h += 1
            continue
        end
        
        rank += 1
        h_carry = h
        
        # rows of U swapped
        temp = U[i, :]
        U[i, :] = U[pivot_row, :]
        U[pivot_row, :] = temp
        
        pivot_val = U[i, h_carry]
        @inbounds @simd for j in i+1:1:sz[1]
            mbp = U[j, h_carry] / pivot_val 
            @inbounds @simd for k = 1:1:sz[2]
                U[j, k] = U[j, k] - mbp * U[i, k]
            end
        end
        
        i += 1
        h += 1
    end
    return U, rank
end

function nullspace_eigen(U::Matrix{<:Number}; tol::Float64=1e-10)
    sz = size(U)    
    row_status = zeros(Int64, sz[1], 2) 
    column_status = zeros(Int64, sz[2])
    
 @inbounds for i in 1:1:sz[1]
    @inbounds for j in 1:1:sz[2]
            # Check against tolerance instead of != 0
            if abs(U[i,j]) > tol 
                row_status[i,1] = 1
                row_status[i,2] = j
                column_status[j] = 1
                break
            end
        end
    end 

    linear_combination_of_vectors = Matrix{ComplexF64}[]   
    
 @inbounds for j in 1:1:sz[2]
        if column_status[j] == 0 
            # Free column found hence Construct basis vector
            basis_vector = zeros(ComplexF64, sz[2],1)
            basis_vector[j] = 1.0 + 0.0im 
            
          @inbounds for i in sz[1]:-1:1
                if row_status[i,1] != 0
                    pivot_column = row_status[i,2]
                    b_comp = 0.0 + 0.0im # RHS is always 0
                    
                    @inbounds @simd for doraemon in pivot_column+1:1:sz[2]
                        b_comp -= U[i,doraemon] * basis_vector[doraemon]
                    end
                    basis_vector[pivot_column] = b_comp / U[i,pivot_column]
                end
            end
            push!(linear_combination_of_vectors, basis_vector)
        end
    end 

    return linear_combination_of_vectors   
end

function eigenvectors(A::Matrix{<:Number}; tol=1e-10)
    n = size(A, 1)
    eigenvals = eigenvalues(A)
    all_eigenvectors = Any[]
    
    for lamda in eigenvals
        #  Construct (A - lamda*I)
        D = ComplexF64.(A)
        for i in 1:n
            D[i, i] -= lamda
        end
        
        # 2. Reduce using custom eigensolver GE
        U, rank = ref_eigen(D, tol=tol)
        
        # 3. Extract the null space basis vectors
        if rank < n
            e_vecs = nullspace_eigen(U, tol=tol)
            append!(all_eigenvectors, e_vecs)
        else
            append!(all_eigenvectors, ["Failed to find null space. Adjust tolerance!"])
        end
    end
    
    return all_eigenvectors
end