function add(u::Vector{<:Number},v::Vector{<:Number})
        ulen = size(u)[1]
        vlen = size(v)[1]
        if ulen != vlen 
            throw(ArgumentError("dimension mismatch"))
        else 
            sum = Vector{promote_type(eltype(u),eltype(v))}(undef, ulen)
           @inbounds @simd for i = 1:1:ulen
                sum[i] = v[i]+u[i]
                           end
        return sum
        end
end

function sub(u::Vector{<:Number},v::Vector{<:Number})
        ulen = size(u)[1]
        vlen = size(v)[1]
        if ulen != vlen 
            throw(ArgumentError("dimension mismatch"))
        else 
            sum = Vector{promote_type(eltype(u),eltype(v))}(undef, ulen)
           @inbounds @simd for i = 1:1:ulen
                sum[i] = u[i]-v[i]
                           end
        return sum
        end
end

function scale(u::Vector{<:Number},sf::Number)
    scaled_vector = Vector{promote_type(eltype(u),eltype(sf))}(undef,size(u)[1])
    @inbounds @simd for i in 1:1:size(u)[1]
        scaled_vector[i] = sf*u[i]
    end
    return scaled_vector
end

function mag(u::Vector{<:Number})
    a = 0
    @inbounds @simd for i in 1:1:size(u)[1]
        a += u[i]^2
    end
    return a^0.5
end

function dot(u::Vector{<:Number},v::Vector{<:Number})
    ulen = size(u)[1]
    vlen = size(v)[1]
    if ulen != vlen
        throw(ArgumentError("Dimension Mismatch"))
    else 
        a = 0
        for i in 1:1:ulen
            a += u[i]*v[i]
        end
    return a
    end
end