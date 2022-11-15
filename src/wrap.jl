function wrap!(box, coords, images=nothing)
    boxlo, boxhi = box[ : , 1], box[ : , 2]
    L = boxhi - boxlo
    if images !== nothing
        images += fld.(coords .- boxlo, L)
    end
    coords = ((coords .- boxlo) .% L) .+ boxlo
    return nothing
end

function unwrap!(box, coords, images, bonds=nothing)  # Not yet implemented
    boxlo, boxhi = box[ : , 1], box[ : , 2]
    L = boxhi - boxlo
    return nothing
end