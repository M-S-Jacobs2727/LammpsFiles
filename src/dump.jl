
struct DumpFrame
    timestep::Int
    natoms::Int
    properties::Vector{String}
    box::Matrix{<:Real}
    atoms::Matrix{<:Real}  # with size (length(properties), natoms) (i.e., each row is a property, each column is an atom)
    idtoindex::Vector{<:Integer}
end

function removeComment(line::AbstractString)
    strip(split(line, "#")[1])
end

# TODO: Improve performance of large dumpfiles using file pointers (like Frank), and
# investigate mmap. -> moving to separate functions, `load_dump` and `read_frame`.
# Alternately, just return an iterator over each snapshot, in order. Doing so will
# require tasks. The main task will interact with the user, and the background task 
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
            return DumpFrame(timestep, natoms, properties, box, atoms, idtoindex)
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
        box = zeros(Float32, 3, 2)
        values = split(removeComment(readline(f)))
        box[1, : ] = [parse(Float32, values[1]), parse(Float32, values[2])]
        values = split(removeComment(readline(f)))
        box[2, : ] = [parse(Float32, values[1]), parse(Float32, values[2])]
        values = split(removeComment(readline(f)))
        box[3, : ] = [parse(Float32, values[1]), parse(Float32, values[2])]
        
        # Atoms
        while !occursin("ITEM: ATOMS", line)
            line = removeComment(readline(f))
        end
        properties = split(line)[3:end]

        atoms = zeros(Float32, length(properties), natoms)  # transposed for performance(?)
        for i=1:natoms
            line = readline(f)
            atoms[:, i] = parse.(Float32, split(line)[1:length(properties)])
        end
    end
    atom_ids = atoms[findfirst(x->x=="id", properties), : ]
    idtoindex = zeros(maximum(atom_ids))
    for (i, a) in enumerate(atom_ids)
        idtoindex[a] = i
    end

    return DumpFrame(timestep, natoms, properties, box, atoms, idtoindex)
end
