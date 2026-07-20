function LU(M :: AbstractMatrix{<:Number})
    U = Float64.(M)
    sz = size(U)
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
        segment = @view U[i+1:sz[1], h]
    @inbounds @simd for meow in 1:1:sz[1]-i
            temp = abs(segment[meow]) 
            if temp > pivot # partial pivoting is occuring here
                pivot = temp
                pivot_row = meow + i
            end
        end
        
        # If the largest element in this column is 0 then skip to the next column
        # Because a number like 0 or near zero leads to large numbers ( when divided by 0 or near zero)
        # which leads to overflow and roundoff errors
        if pivot == 0
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
        end
        
        i += 1
        h += 1
    end
    return P, L, U, rank
end

function Householder_matrix(v::Vector{<:Number}) 
    # it takes a vector and returns a corresponding householder matrix 
    sz = size(v)
    mag2 = 0
    I = zeros(Float64,sz[1],sz[1])
  @inbounds @simd  for i in 1:1:sz[1]  # construction of Identity matrix
        mag2 += v[i]^2
        I[i,i]=1
    end
    if v[1]>=0
        v[1] = v[1]+mag2^0.5
    else
        v[1] = v[1]-mag2^0.5  
    end   
    vTv = (transpose(v)*v)
    vvT = v*transpose(v)
    H = I - 2*(vvT/vTv) # construction of householder matrix

    return H
end
function QR(M_raw :: AbstractMatrix{<:Number})
    R = Float64.(M_raw)
    sz = size(R)
    pikachu = min(sz[1]-1,sz[2])  # max amount of rank possible
    #= diagonal for a non square matrix = minimum of rows and column, and thats exactly whats done here
        although we have modified it a bit for the tall matrices (having more rows then columns)
        now diagonal = minimum of (rows-1) and column, otherwise few columns would not be transformed accordingly
    =#
    Q = zeros(Float64,sz[1],sz[1])
 @inbounds @simd   for i in 1:1:sz[1]  # construction of Identity matrix { Q =I}
        Q[i,i]=1
    end
    
 @inbounds @simd   for i in 1:1:pikachu
        column_vec = @view R[i:sz[1],i] # selcting the column vector which is about to get boom-bam budup-budup booooooooooom
                                        # householder matrix will be constructed according to this vector
        #------------------------------------------------------------------------------------------------------
        hehe = Householder_matrix(copy(column_vec)) 
        R[i:sz[1],i:sz[2]] = hehe*R[i:sz[1],i:sz[2]] 
        #=householder matrix is applied to the effective matrix
        a householder matrix transforms the column vector (which is already selected above)
        hence the effective application of householder matrix is from the ith row to last row and ith column to last column
        {as the matrix before ith row and column is already converted into an upper triangular matrix)
        =#
       
        Q[:, i:sz[1]] = Q[:, i:sz[1]] * hehe # updating the Q matrix => Q = H1*H2*H3.......Hn
        
        #---------------------------------------------------------------------------------------------
        R[i+1:sz[1],i] .= 0.0 # without this also the algo will work, but it will have floating point errors
        # to make the elements below diagonal as zero, broadcasting{dot operator} "." is done
        # it is not necessary but it makes the matrix pleasing to eyes
    end
    return Q,R
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
                    b_comp = 0
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