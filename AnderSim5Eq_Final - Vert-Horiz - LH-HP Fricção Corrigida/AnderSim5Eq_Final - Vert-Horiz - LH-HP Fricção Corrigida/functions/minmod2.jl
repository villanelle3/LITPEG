function minmod2(v::Matrix{Float64})
    m, n = size(v)
    mfunc = zeros(1, n)
    s = sum(sign.(v), dims=1) ./ m
    
    @inbounds for j in 1:n
        if abs(s[j]) == 1
            mfunc[j] = s[j] * minimum(abs.(@view v[:, j]))
        end
    end
    
    return mfunc
end
