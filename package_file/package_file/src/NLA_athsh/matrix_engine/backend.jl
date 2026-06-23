function matmult( x:: Matrix{<: Number} , y:: Matrix{<: Number})   
        x_dim =  size(x)
        y_dim =  size(y)
        if x_dim[2] != y_dim[1] 
            throw(ArgumentError("Invalid input detectet/dimensions not matching"))
        else 
            # xy = Matrix{Number}(undef,x_dim[1],y_dim[2])
             xy = zeros(Number,x_dim[1],y_dim[2])
            @inbounds @simd for i in 1:1:y_dim[2] 
           @inbounds @simd     for j in 1:1:x_dim[2] 
               @inbounds @simd   for k in 1:1:x_dim[1] 
                        xy[k,i] += x[k,j]*y[j,i]
                    
                                    end
                                end
                            end
            return xy
        end
 
    end

function add(x:: Matrix{<:Number},y:: Matrix{<: Number})
    x_dim = size(x)
    y_dim = size(y)
    if x_dim != y_dim
        throw(ArgumentError("Dimension Mismatch"))
    else 
        xy = zeros(Number,x_dim[1],x_dim[2])
      @inbounds @simd  for i = 1:1:x_dim[2] # i is  column pointer
         @inbounds @simd   for j = 1:1:x_dim[1] # j is row pointer
                xy[j,i] = x[j,i]+y[j,i]
            end
        end       
    end
    return xy
end

function sub(x:: Matrix{<:Number},y:: Matrix{<: Number})
    x_dim = size(x)
    y_dim = size(y)
    if x_dim != y_dim
        throw(ArgumentError("Dimension Mismatch"))
    else 
        xy = zeros(Number,x_dim[1],x_dim[2])
      @inbounds @simd  for i = 1:1:x_dim[2] # i is  column pointer
         @inbounds @simd   for j = 1:1:x_dim[1] # j is row pointer
                xy[j,i] = x[j,i] - y[j,i]
            end
        end       
    end
    return xy
end

function transpose_(x::Matrix{<:Number})
    dims = size(x)
    xT =zeros(Number, dims[2],dims[1])
    @inbounds @simd for i = 1:1:dims[1] 
       @inbounds @simd for j in 1:1:dims[2]  
                              xT[j,i] = x[i,j]
        end
    end
return xT
end

function trace(M_raw ::Matrix{<:Number})
    M= M_raw
    sz = size(M)

    if sz[2] != sz[1]
        return throw(ArgumentError("Non-Square Matrices dont have Trace"))
    else
        trace = 0
        for i in 1:sz[1]
            trace += M[i, i]
        end
    end
    return trace
end




