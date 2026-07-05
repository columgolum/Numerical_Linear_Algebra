function ref(M :: Matrix{<:Number})
    U = Float64.(M)
    sz =size(U)
    i =1
    rank = 0
 while i<= sz[1] # Selects/locks the rows(selection of pivot row), The pivot is selected in this loop onlu
            pivot = 0
            h_carry = 0
            h =1
  while h<=sz[2]
                if U[i,h] != 0
                    pivot = U[i,h]
                    h_carry = h
                    rank+=1
                    break
                else
                     meow = i+1
                    if meow <= sz[1]
            while meow<= sz[1]
                            if  U[meow,h] != 0
                            
                                temp = U[meow, :]
                                U[meow, :] = U[i, :]
                                U[i, :] = temp
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
                end
            end
            i+=1
    end
    return U,rank,sz[2]-rank
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

function solveLSE(D::Matrix{<:Number})
    M ,rank,nullity = rref(D)
    sz = size(M)
    for i in 1:sz[1]
        count = 0 
        for j in 1:(sz[2]-1)
            if M[i, j] != 0
                count += 1
                break # if pivot found then the row will give a valid answer, hence break and start from next row
            end
        end
        
        # in A|b : b = M[i, sz[2]]
        # if b is non zero, but the matrix A has zero rows then it means no solution
        if count == 0 && M[i, sz[2]] != 0
            return ("No solution exists",) # returns a tuple which is resolved accordingly at the python junction
        end
    end
    if rank < (sz[2]-1)
        particular_vector ,linear_combination_of_vectors = solutions(M,rank)
        return "infinite soln" , particular_vector ,linear_combination_of_vectors
    else
        return "Unique solution",solutions(M,rank) 
    end
end
#-------------------------------------------------------------------------------------------------------------------------------------------

function solutions(D::Matrix{<:Number},r::Number) # D = rref(A)
    M_raw = D
    rank = r
    sz = size(M_raw)    
    M = Float64.(M_raw[:,1:sz[2]-1])   
    
    column_status = zeros(Int64, sz[2]-1, 2) # pivot column = 1 , free column = 0
    # first Column  = keeps the check on pivot columns
    # second column = keeps the track in which that pivot resides
    
  @inbounds for i in 1:1:sz[1]
  @inbounds for j in 1:1:sz[2]-1
            if M[i,j] != 0
                column_status[j,1] = 1
                column_status[j,2] = i
                break
            end
        end
    end #= By the end of this code block, our column_status Matrix(or array) should have the information of pivots in which column(in binary)
    and in which row =#

    
    particular_vector = zeros(Float64,sz[2]-1,1) 
 @inbounds @simd   for p_col in 1:1:sz[2]-1 # traverse the columns
        if column_status[p_col, 1] == 1     # if a pivot is found
            row = column_status[p_col, 2]   # then look for row in which that resides 
            particular_vector[p_col] = M_raw[row, sz[2]]  # then the particular vector will be the last element of that row

            # all free variables are set to be 0
        end
    end # particular vector is done after this code block

    
    linear_combination_of_vectors = Matrix{Float64}[]   
@inbounds @simd    for j in 1:1:sz[2]-1
            if column_status[j,1] == 0 #= We only care about free columns where column 1 is 0 (i.e column_status[j,1] = 0)
                                         Because the pivot variabels of basis vector are (-1)*M[:,j] where j is a free column
                                        =#
            
            basis_vector = zeros(Float64,sz[2]-1,1)
            basis_vector[j] = 1 # the variable of free column(free variable is set to 1 => this is linear algebra theory)
            
            # Now, look at all the pivot columns to fill in the rest of this vector
        @inbounds @simd    for p_col in 1:1:sz[2]-1 # this loop traverses the columns 
                if column_status[p_col,1] == 1 
                    row = column_status[p_col, 2] #  note down the rows in which the pivot resides
                    basis_vector[p_col] = -M[row, j] #=  assign the element from the free column(after flipping the sign) 
                                                       to the pivot column element of the basis vector =#
                end
            end
            
            push!(linear_combination_of_vectors, basis_vector)
        end
    end
    if rank == sz[2]-1
         return particular_vector
    else
         return (particular_vector ,linear_combination_of_vectors)
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
function det(M_raw ::Matrix{<:Number},)
    M,rank, nullity = ref(M_raw)
    sz = size(M)

    if sz[2] != sz[1]
        return throw(ArgumentError("Non-Square Matrices dont have determinant"))
    elseif  rank<sz[2]
        return 0
    else
        determinant = 1
        for i in 1:sz[1]
            determinant *= M[i, i]
        end
    end
    return determinant
end
