function ref(M :: Matrix{<:Number})
    U = Float64.(M)
    sz = size(U)
    i = 1
    rank = 0
    h = 1 
    while i <= sz[1] && h <= sz[2] # traversing both rows and columns
        pivot_row = i
        pivot = abs(U[i, h])
        
        # Scan ALL rows below row i in the current column h
        segment = @view U[i+1:sz[1], h]
      @inbounds @simd  for meow in 1:1:sz[1]-i
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
    return U, rank,sz[2]-rank
end
#------------------------------------------------------------------------------------------------------------------------------

function rref(U::Matrix{<:Number})
     D,rank,nullity = ref(U)
    
    sz = size(D)
 @inbounds  for i = sz[1]:-1:1 # goes from last row to the first row
        pivot = 0
        h_carry = 0
   @inbounds   for h in 1:1:sz[2] # traverses the selected row
            if D[i,h] != 0 # if a pivot is found then 
                h_carry = h
                break
            end
        end
        if h_carry != 0 # if h_carry is there then a pivot is found, hence only then we will proceed to do backward( OR upward) elimination 
            pivot = D[i,h_carry] 
      @inbounds @simd      for meow in 1:1:sz[2]
                D[i,meow] = D[i,meow]/pivot
            end
      @inbounds @simd for j in i-1:-1:1 # selects/ locks the rows above the pivot row( selected by i pointer) 
                    mbp = D[j,h_carry]
               @inbounds @simd     for k in 1:1:sz[2]
                        D[j,k] = D[j,k]-mbp*D[i,k]
                    end
            end 
        end
    end
    return D,rank,nullity
end




#--------------------------------------------------------------------------------------------------------------------------------------------

function solveLSE(D::Matrix{<:Number}, mode::String)
    sz = size(D)
    
    # Take out the matrix A and vector b from the augmented matrix D
    nov = sz[2] - 1 # nov = number of variables
    A = @view D[:, 1:nov]
    b = copy(D[:, sz[2]])

    if mode == "GE"
        M, rank, nullity = ref(D)
        for i in 1:sz[1]
            count = 0 
            for j in 1:1:nov
                if M[i, j] != 0
                    count += 1
                    break 
                end
            end
            
            # Check if in A|b the last column is non-zero while all other columns are zero
            if count == 0 && M[i, sz[2]] != 0
                return ("No solution exists",) 
            end
        end
        
        if rank < nov
            particular_vector, linear_combination_of_vectors = solutions(M, rank)
            return "infinite soln", particular_vector, linear_combination_of_vectors
        else
            return "Unique solution", solutions(M, rank) 
        end
        
    elseif mode == "LU"
        P, L, U, rank = LU(A)
        sz_of_L = size(L)
        b_new = P * b 
        y = zeros(Float64, sz_of_L[2], 1)
        
        @inbounds for i in 1:1:sz_of_L[1]
            b_comp = b_new[i] 
            @inbounds @simd  for doraemon in i-1:-1:1
                b_comp -= L[i, doraemon] * y[doraemon]
            end
            y[i] = b_comp
        end  
        
        for i in 1:sz[1]
            count = 0 
            for j in 1:1:nov
                if U[i, j] != 0
                    count += 1
                    break 
                end
            end
            
            if count == 0 && y[i] != 0
                return ("No solution exists",) 
            end
        end
        # ------------------------------
        
        if rank < nov
            particular_vector, linear_combination_of_vectors = solutions(U, y, rank) 
            return "infinite soln", particular_vector, linear_combination_of_vectors
        else
            return "Unique solution", solutions(U, y, rank) 
        end
        
    elseif mode == "QR"
        Q, R = QR(A)
        
        b_new = transpose(Q) * b
        
        for i in 1:sz[1]
            count = 0 
            for j in 1:1:nov
                if R[i, j] != 0
                    count += 1
                    break 
                end
            end
            
            if count == 0 && b_new[i] != 0
                return ("No solution exists",) 
            end
        end
        # ------------------------------
        rank = 0
        for i in 1:1:sz[1]
            row_vec = @view R[i, 1:nov]
            for j in 1:1:nov
                if row_vec[j] == 0 
                    nothing
                else
                    rank += 1
                    break
                end
            end
        end

        if rank < nov 
            particular_vector, linear_combination_of_vectors = solutions(R, b_new, rank) 
            return "infinite soln", particular_vector, linear_combination_of_vectors
        else
            return "Unique solution", solutions(R, b_new, rank) 
        end
        
    else
        throw(ArgumentError("Supported modes are = GE, LU, QR"))
    end
end
#-------------------------------------------------------------------------------------------------------------------------------------------

function solutions(D::Matrix{<:Number},r::Number) # D = ref(A) : D is augmented matrix
    U = D
    rank = r
    sz = size(U)    

    particular_vector = zeros(Float64,sz[2]-1,1)
    row_status = zeros(Int64,sz[1],2) 
    # rows of row_status = number of rows
    # first column tells wheather a row is pivot row(1) or free row(0) 
    # second column tells in which column that pivot is 
    column_status = zeros(Int64,sz[2]-1,1)
    # columns of column_status = number of columns in A which is A|b = U
    # first row tells wheather a column is pivot column(1) or free column(0)
    
 @inbounds  for i in 1:1:sz[1]
    @inbounds     for j in 1:1:sz[2]-1
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
            b_comp = U[i,sz[2]]
          @inbounds @simd  for doraemon in pivot_column+1:1:sz[2]-1
                if column_status[doraemon] == 0
                    nothing
                else
                    b_comp -= U[i,doraemon]*particular_vector[doraemon]
                end
            end
            particular_vector[pivot_column]= b_comp/U[i,pivot_column]
        end
    end  # this code block is used to find the particular solution of the system of linear equations, and store it in particular_vector

    if rank == sz[2]-1
        return particular_vector
    else 

       linear_combination_of_vectors = Matrix{Float64}[]   
@inbounds for j in 1:1:sz[2]-1
            if column_status[j] == 0 #= We only care about free columns where column_status[j] == 0
                                         Because the pivot variabels of basis vector are (-1)*M[:,j] where j is a free column
                                        =#
            
            basis_vector = zeros(Float64,sz[2]-1,1)
            basis_vector[j] = 1 # the variable of free column(free variable is set to 1 => this is linear algebra theory)
          @inbounds for i in sz[1]:-1:1
                if row_status[i,1] ==0
                    nothing
                else
                    pivot_column = row_status[i,2]
                    b_comp = 0
                    @inbounds @simd  for doraemon in pivot_column+1:1:sz[2]-1
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
#------------------------------------------------------------------------------------------------------------------------------------------------
function inverse(M::Matrix{<:Number})
    sz = size(M)
     U = Float64.(M)
    if sz[1] != sz[2]
        return throw(ArgumentError("Non Square Matrices dont have inverse"))
    #= elseif ref(M)[2] <sz[2]
        return "Inverse do NOT exist" =#
    elseif ref(M)[2] < sz[2]
        return "Inverse Do Not Exist"
    else
       I_mat =  zeros(Float64,sz[1],sz[2])
        for i in 1:1:sz[1]
            I_mat[i,i] = 1
        end
          i =1
 while i<= sz[1] # Selects/locks the rows(selection of pivot row), The pivot is selected in this loop onlu
            pivot = 0
            h_carry = 0
            h =1
  while h<=sz[2]
                if U[i,h] != 0
                    pivot = U[i,h]
                    h_carry = h
                    break
                else
                    meow = i+1
                    
                    if meow <= sz[1]
            while meow<= sz[1]
                            if  U[meow,h] != 0
                            
                                temp1 = U[meow, :]
                                U[meow, :] = U[i, :]
                                U[i, :] = temp1

                                temp2 = I_mat[meow, :]
                                I_mat[meow, :] = I_mat[i, :]
                                I_mat[i, :] = temp2
                                break                                   
                            else
                                if meow != sz[1]
                                    meow+=1
                                else # meow == sz[1]
                                    h+=1
                                    break
                                end
                            end
                        end
                    else
                        h+=1
                    end
                end
            end  
            if pivot == 0
                i+=1
                continue
            end   
@inbounds @simd  for j in i+1:1:sz[1] # its selects rows on which row operation is about to occur, multiplier is also selected in this loop
                mbp = U[j,h_carry]/pivot # mbp = multiplier by pivot i.e multiplier/pivot
@inbounds @simd for k = 1:1:sz[2]
                    U[j,k] = U[j,k]-mbp*U[i,k]
                     I_mat[j,k] = I_mat[j,k]-mbp*I_mat[i,k]
                end
            end
            i+=1
    end
    for i = sz[1]:-1:1 # goes from last row to the first row
        pivot = 0
        h_carry = 0
        for h in 1:1:sz[2] # traverses the selected row
            if U[i,h] != 0 # if a pivot is found then 
                h_carry = h
                break
            end
        end
        if h_carry != 0 # if h_carry is there then a pivot is found, hence only then we will proceed to do backward( OR upward) elimination 
            pivot = U[i,h_carry] 
            for meow in 1:1:sz[2]
                U[i,meow] = U[i,meow]/pivot
                I_mat[i,meow] = I_mat[i,meow]/pivot
            end
            for j in i-1:-1:1 # selects/ locks the rows above the pivot row( selected by i pointer) 
                    mbp = U[j,h_carry]
                    for k in 1:1:sz[2]
                        U[j,k] = U[j,k]-mbp*U[i,k]
                        I_mat[j,k] = I_mat[j,k]-mbp*I_mat[i,k]
                    end
            end 
        end
    end
    return I_mat
              
end  
end
#----------------------------------------------------------------------------------------------------------------------------------------------------
function det(M::Matrix{<:Number})
    U = Float64.(M)
    sz = size(U)
    if sz[2] != sz[1]
        return throw(ArgumentError("Non-Square Matrices dont have determinant"))
    end
    i = 1
    swaps = 0
    h = 1
    rank = 0
    while i <= sz[1] && h <= sz[2] # traversing both rows and columns
        pivot_row = i
        pivot = abs(U[i, h])
        
        # Scan ALL rows below row i in the current column h
        segment = @view U[i+1:sz[1], h]
      @inbounds @simd  for meow in 1:1:sz[1]-i
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
        if pivot_row != i
            temp = U[i, :]
            U[i, :] = U[pivot_row, :]
            U[pivot_row, :] = temp
            swaps += 1
        end

        
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
    
    if  rank<sz[2]
        return 0
    else
        determinant = 1
        for i in 1:sz[1]
            determinant *= U[i, i]
        end
    end
    return determinant*(-1)^swaps
end
