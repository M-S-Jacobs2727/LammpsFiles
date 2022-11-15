using Logging

function removeComment(line::AbstractString)
    strip(split(line, "#")[1])
end

# TODO: Improve performance of large dumpfiles via tasks and an iterator
# The main task will interact with the user, and the background task 
# will be reading and holding the place in the file.
# 
# TODO: Implement dump_style yaml.
"""
Read a single frame of data from a dump file.
Currently only supports dump styles `custom` and `atom`.
"""
function readDump(source::AbstractString)
    timestep = 0
    natoms = 0
    properties = Vector{String}(undef, 1)
    box = zeros(3, 2)
    tilt = zeros(3)
    atoms = zeros(1, 1)
    idtoindex = zeros(1)
    open(source, "r") do f
        line = ""

        # Timestep, or truncated file
        while !eof(f)
            line = readline(f)
            if removeComment(line) == "ITEM: TIMESTEP"
                break
            end
        end
        if eof(f)
            @warn "Reached end of file before finding anything!"
            return (
                timestep=timestep,
                natoms=natoms,
                properties=properties,
                box=box,
                tilt=tilt,
                atoms=atoms,
                idtoindex=idtoindex
            )
        end
        timestep = parse(Int, removeComment(readline(f)))

        while removeComment(line) != "ITEM: NUMBER OF ATOMS"
            line = readline(f)
        end
        natoms = parse(Int, removeComment(readline(f)))
        
        # Box
        while !occursin("ITEM: BOX BOUNDS", removeComment(line))
            line = readline(f)
        end
        box = zeros(Float64, 3, 2)
        tilt = zeros(Float64, 3)
        hastilt = occursin("xy", line)
        
        values = split(removeComment(readline(f)))
        box[1, : ] = [parse(Float64, values[1]), parse(Float64, values[2])]
        if hastilt tilt[1] = parse(Float64, values[3]) end
        
        values = split(removeComment(readline(f)))
        box[2, : ] = [parse(Float64, values[1]), parse(Float64, values[2])]
        if hastilt tilt[2] = parse(Float64, values[3]) end
        
        values = split(removeComment(readline(f)))
        box[3, : ] = [parse(Float64, values[1]), parse(Float64, values[2])]
        if hastilt tilt[3] = parse(Float64, values[3]) end
        
        # Atoms
        while !occursin("ITEM: ATOMS", line)
            line = removeComment(readline(f))
        end
        properties = split(line)[3:end]

        atoms = zeros(Float64, length(properties), natoms)
        for i in axes(atoms, 2)
            line = readline(f)
            atoms[:, i] = parse.(Float64, split(line)[1:length(properties)])
        end
    end

    atom_ids = atoms[findfirst(x->x=="id", properties), : ]
    idtoindex = zeros(maximum(atom_ids))
    for (i, a) in enumerate(atom_ids)
        idtoindex[a] = i
    end

    return (
        timestep=timestep,
        natoms=natoms,
        properties=properties,
        box=box,
        tilt=tilt,
        atoms=atoms,
        idtoindex=idtoindex
    )
end
