"""
    wrap!(coords::Matrix{<:Real}, box::Matrix{<:Real})

Wrap molecular coordinates into a box with bounds
`box = [xlo xhi; ylo yhi; zlo zhi]`. Returns the resulting image values,
a Matrix with the same size as `coords`.
"""
function wrap!(coords::Matrix{<:Real}, box::Matrix{<:Real})
    boxlo, boxhi = box[ : , 1], box[ : , 2]
    L = boxhi - boxlo
    images = fld.(coords .- boxlo, L)
    coords = ((coords .- boxlo) .% L) .+ boxlo
    return images
end

"""
    unwrap!(coords::Matrix{<:Real}, box::Matrix{<:Real}, images::Matrix{<:Integer})

Unwrap molecular coordinates from a box with bounds
`box = [xlo xhi; ylo yhi; zlo zhi]` and image values given by `images`,
a Matrix with the same size as `coords`.
"""
function unwrap!(coords::Matrix{<:Real}, box::Matrix{<:Real}, images::Matrix{<:Integer})
    L = box[ : , 2] - box[ : , 1]
    coords += L .* images
    return nothing
end