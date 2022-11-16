using Logging

function removeComment(line::AbstractString)
    strip(split(line, "#")[1])
end

function readDataSection(f::IOStream, data::Matrix{Any})
    datatype = eltype(data)
    for i in axes(data, 2) #eachindex(data[1, : ])
        line = readline(f)
        data[ : , i] = parse.(datatype, split(removeComment(line)))
    end
    return data
end

"""
Unsupported atom_styles: template, hybrid, spin, dielectric, dipole
"""
function readData(source::AbstractString, atom_style::AbstractString="full")
    if atom_style == "atomic"
        base_num_cols = 5
    elseif (
        atom_style == "molecular"
        || atom_style == "angle"
        || atom_style == "bond"
        || atom_style == "charge"
        || atom_style == "dpd"
        || atom_style == "mdpd"
        )
        base_num_cols = 6
    elseif (
        atom_style == "full"
        || atom_style == "body"
        || atom_style == "edpd"
        || atom_style == "ellipsoid"
        || atom_style == "peri"
        || atom_style == "sphere"
        )
        base_num_cols = 7
    elseif (
        atom_style == "electron"
        || atom_style == "line"
        || atom_style == "sph"
        || atom_style == "bpm/sphere"
        || atom_style == "template"
        || atom_style == "tri"
        )
        base_num_cols = 8
    elseif atom_style == "mesont" || atom_style == "wavepacket"
        base_num_cols = 11
    elseif atom_style == "smd"
        base_num_cols = 13
    else
        throw(ArgumentError("invalid atom_style `$(atom_style)`"))
    end
    
    return_dict = Dict{String, <:Array{<:Real}}()
    open(source, "r") do f
        readline(f)  # skip first line
        line = readline(f)  # next line blank

        # Header values
        natoms = natom_types = nbonds = nbond_types = nangles = nangle_types = \ 
            ndihedrals = ndihedral_types = nimpropers = nimproper_types = 0

        image_flag = false

        # Parse header, read until non-blank line without header keyword
        while !eof(f)
            line = readline(f)
            line = removeComment(line)
            if line == ""
                continue
            end
            
            if occursin("atoms", line)
                natoms = parse(Int, split(line)[1])
            elseif occursin("bonds", line)
                nbonds = parse(Int, split(line)[1])
            elseif occursin("angles", line)
                nangles = parse(Int, split(line)[1])
            elseif occursin("dihedrals", line)
                ndihedrals = parse(Int, split(line)[1])
            elseif occursin("impropers", line)
                nimpropers = parse(Int, split(line)[1])
            elseif occursin("atom types", line)
                natom_types = parse(Int, split(line)[1])
            elseif occursin("bond types", line)
                nbond_types = parse(Int, split(line)[1])
            elseif occursin("angle types", line)
                nangle_types = parse(Int, split(line)[1])
            elseif occursin("dihedral types", line)
                ndihedral_types = parse(Int, split(line)[1])
            elseif occursin("improper types", line)
                nimproper_types = parse(Int, split(line)[1])
            elseif occursin("xlo xhi", line)
                box = zeros(3, 2)
                values = split(line)
                box[1, : ] = [parse(Float64, values[1]), parse(Float64, values[2])]
                values = split(removeComment(readline(f)))
                box[2, : ] = [parse(Float64, values[1]), parse(Float64, values[2])]
                values = split(removeComment(readline(f)))
                box[3, : ] = [parse(Float64, values[1]), parse(Float64, values[2])]
                return_dict["box"] = box
            elseif occursin("xy xz yz", line)
                values = split(line)
                return_dict["tilt"] = [parse(Float64, values[1]), parse(Float64, values[2]), parse(Float64, values[3])]
            elseif occursin("extra bond per atom", line)
                continue
            elseif occursin("extra angle per atom", line)
                continue
            elseif occursin("extra dihedral per atom", line)
                continue
            elseif occursin("extra improper per atom", line)
                continue
            elseif occursin("extra special per atom", line)
                continue
            elseif occursin("ellipsoids", line)
                continue
            elseif occursin("lines", line)
                continue
            elseif occursin("triangles", line)
                continue
            elseif occursin("bodies", line)
                continue
            else
                break
            end
        end

        # Body sections supported: 
        #       Atoms, Velocities, Bonds, Angles, Dihedrals, Impropers
        # Body sections not supported (per-atom-type instead of per-atom):
        #       Masses, Pair Coeffs, PairIJ Coeffs, Bond Coeffs, Angle Coeffs,
        #       Dihedral Coeffs, Improper Coeffs
        while !eof(f)
            if removeComment(line) == ""
                line = readline(f)
                continue
            end
            if occursin("Atoms", line)
                m = match(r"^Atoms.*(# +([a-z]+))", line)
                if m !== nothing
                    atom_style_hint = m.captures[2]
                    if atom_style != atom_style_hint
                        @warn "Provided atom_style $(atom_style) does not match atom_style hint \
                            found in data file, $(atom_style_hint). Using provided atom_style $(atom_style)"
                    end
                end
                readline(f)  # Skip next line

                line = readline(f)
                first_row_data = parse.(Float64, split(removeComment(line)))
                
                num_cols = length(first_row_data)
                if num_cols != base_num_cols && num_cols != base_num_cols + 3
                    throw(DimensionMismatch("invalid number of columns in Atoms section of data file: \
                                         $(num_cols). Expected $(base_num_cols) or $(base_num_cols+3)"))
                end
                if num_cols == base_num_cols + 3
                    image_flag = true
                end

                tmp_data = zeros(num_cols, return_dict["natoms"])
                tmp_data[ : , 1] = first_row_data
                tmp_data[ : , 2:end] = readDataSection(f, tmp_data[ : , 2:end])

                indices = sortperm(tmp_data[1, : ])
                tmp_data = tmp_data[ : , indices]
                
                return_dict["atomIDs"] = Int.(tmp_data[1, : ])
                if atom_style == "atomic"
                    return_dict["atom_types"] = Int.(tmp_data[2, : ])
                    return_dict["coords"] = tmp_data[3:5, : ]'
                elseif atom_style == "charge"
                    return_dict["atom_types"] = Int.(tmp_data[2, : ])
                    return_dict["charges"] = tmp_data[3, : ]
                    return_dict["coords"] = tmp_data[4:6, : ]'
                elseif (
                    atom_style == "molecular"
                    || atom_style == "bond"
                    || atom_style == "angle"
                )
                    return_dict["molecules"] = Int.(tmp_data[2, : ])
                    return_dict["atom_types"] = Int.(tmp_data[3, : ])
                    return_dict["coords"] = tmp_data[4:6, : ]'
                elseif atom_style == "full"
                    return_dict["molecules"] = Int.(tmp_data[2, : ])
                    return_dict["atom_types"] = Int.(tmp_data[3, : ])
                    return_dict["charges"] = tmp_data[4, : ]
                    return_dict["coords"] = tmp_data[5:7, : ]'
                elseif (
                    atom_style == "line"
                    || atom_style == "mesont"
                    || atom_style == "bpm/sphere"
                    || atom_style == "tri"
                )
                    return_dict["molecules"] = Int.(tmp_data[2, : ])
                    return_dict["atom_types"] = Int.(tmp_data[3, : ])
                    return_dict["the_rest"] = image_flag ? tmp_data[4:end-6, : ]' : tmp_data[4:end-3, : ]'
                    return_dict["coords"] = image_flag ? tmp_data[end-5:end-3, : ]' : tmp_data[end-2:end, : ]'
                else
                    return_dict["atom_types"] = Int.(tmp_data[2, : ])
                    return_dict["the_rest"] = image_flag ? tmp_data[4:end-6, : ]' : tmp_data[4:end-3, : ]'
                    return_dict["coords"] = image_flag ? tmp_data[end-5:end-3, : ]' : tmp_data[end-2:end, : ]'
                end
                if image_flag
                    return_dict["images"] = Int.(tmp_data[end-2:end, : ]')
                end
                idtoindex = zeros(maximum(atomIDs))
                for (i, a) in enumerate(atomIDs)
                    idtoindex[a] = i
                end
                return_dict["idtoindex"] = idtoindex

            elseif occursin("Velocities", line)
                readline(f)  # Blank
                data = zeros(4, natoms)
                data = readDataSection(f, data)
                # Try applying indices from Atoms section to sort velocities by ID
                # Check if sorted. If not, sort normally. Will usually save time.
                issorted(data[1, indices]) || (indices = sortperm(data[1, : ]))
                return_dict["velocities"] = data[2:4, indices]'
            elseif occursin("Masses", line)  # Not yet implemented
                readline(f)  # Blank
                for i in 1:natom_types; readline(f); end
            elseif nbonds > 0 && occursin("Bonds", line)
                readline(f)  # Blank
                data = zeros(Int, 4, nbonds)
                data = readDataSection(f, data)
                return_dict["bond_types"] = data[2, : ]
                return_dict["bonds"] = data[3:4, : ]'
            elseif nangles > 0 && occursin("Angles", line)
                readline(f)  # Blank
                data = zeros(Int, 5, nangles)
                data = readDataSection(f, data)
                return_dict["angle_types"] = data[2, : ]
                return_dict["angles"] = data[3:5, : ]'
            elseif ndihedrals > 0 && occursin("Dihedrals", line)
                readline(f)  # Blank
                data = zeros(Int, 6, ndihedrals)
                data = readDataSection(f, data)
                return_dict["dihedral_types"] = data[2, : ]
                return_dict["dihedrals"] = data[3:6, : ]'
            elseif nimpropers > 0 && occursin("Impropers", line)
                readline(f)  # Blank
                data = zeros(Int, 7, nimpropers)
                data = readDataSection(f, data)
                return_dict["improper_types"] = data[2, : ]
                return_dict["impropers"] = data[3:6, : ]'
            elseif occursin("Pair Coeffs", line)  # Not yet implemented
                readline(f)  # Blank
                for i in 1:natom_types; readline(f); end
            elseif occursin("PairIJ Coeffs", line)  # Not yet implemented
                readline(f)  # Blank
                for i in 1:(natom_types*(natom_types-1)); readline(f); end
            elseif nbonds > 0 && occursin("Bonds Coeffs", line)  # Not yet implemented
                readline(f)  # Blank
                for i in 1:nbond_types; readline(f); end
            elseif nangles > 0 && occursin("Angle Coeffs", line)  # Not yet implemented
                readline(f)  # Blank
                for i in 1:nangle_types; readline(f); end
            elseif ndihedrals > 0 && occursin("Dihedral Coeffs", line)  # Not yet implemented
                readline(f)  # Blank
                for i in 1:ndihedral_types; readline(f); end
            elseif nimpropers > 0 && occursin("Improper Coeffs", line)  # Not yet implemented
                readline(f)  # Blank
                for i in 1:nimproper_types; readline(f); end
            end
            line = readline(f)
        end
    end
    return_dict
end

"""
Currently, very few checks are performed. Should account for atom style, natoms is 
consistent, etc.
"""
function verifySnapshot(snapshot::Dict{String, <:Array{<:Real}}, atom_style::AbstractString)
    atomIDs = snapshot["atomIDs"]
    ndims(atomIDs) != 1 && throw(DimensionMismatch(
        "atomIDs must have 1 dimenstion, found $(ndims(atomIDs))"
        ))
    natoms = length(atomIDs)
    
    atom_types = snapshot["atom_types"]
    ndims(atom_types) != 1 && throw(DimensionMismatch(
        "atom_types must have 1 dimenstion, found $(ndims(atom_types))"
        ))
    length(atom_types) != natoms && throw(DimensionMismatch(
        "atom_types must have the same length as atomIDs,"
        " found $(length(atom_types)) and $natoms, respectively"
    ))
    natom_types = maximum(atom_types)
        
    coords = snapshot["coords"]
    ndims(coords) != 2 && throw(DimensionMismatch(
        "coords must have 1 dimenstion, found $(ndims(coords))"
    ))
    size(coords, 1) != natoms && throw(DimensionMismatch(
        "coords must have the same dim-1 length as atomIDs,"
        " found $(size(coords, 1)) and $natoms, respectively"
    ))
    size(coords, 2) != 3 && throw(DimensionMismatch(
        "coords must have 3 columns, found $(size(coords, 2))"
    ))

end

function writeData(datafile::AbstractString, snapshot::Dict{String, <:Array{<:Real}}, atom_style::AbstractString)
    if !(atom_style in ["atomic", "molecular", "angle", "bond", "charge", "full"])
        throw(ArgumentError(
            "Invalid atom_style '$atom_style'. Must be one of"
            " 'atomic', 'molecular', 'angle', 'bond', 'charge', or 'full'."
        ))
    end

    verifySnapshot(snapshot, atom_style)

    hasbonds = haskey(snapshot, "bonds")  # implies hasmolecules
    hasangles = haskey(snapshot, "angles")
    hasdihedrals = haskey(snapshot, "dihedrals")
    hasimpropers = haskey(snapshot, "impropers")
    hastilt = haskey(snapshot, "tilt")
    hasimages = haskey(snapshot, "images")
    hascharges = haskey(snapshot, "charges")
    hasvelocities = haskey(snapshot, "velocities")

    natoms = length(snapshot["atom_types"])
    natom_types = maximum(snapshot["atom_types"])
    if hasbonds
        nbonds = length(snapshot["bond_types"])
        nbond_types = maximum(snapshot["bond_types"])
    end
    if hasangles
        nangles = length(snapshot["angle_types"])
        nangle_types = maximum(snapshot["angle_types"])
    end
    if hasdihedrals
        ndihedrals = length(snapshot["dihedral_types"])
        ndihedral_types = maximum(snapshot["dihedral_types"])
    end
    if hasimpropers
        nimpropers = length(snapshot["improper_types"])
        nimproper_types = maximum(snapshot["improper_types"])
    end

    open(datafile, "w") do f
        write(f, "# LAMMPS data file written by Michael Jacobs' LammpsFiles code\n\n")
        write(f, "$natoms atoms\n")
        if hasbonds write(f, "$nbonds bonds") end
        if hasangles write(f, "$nangles angles") end
        if hasdihedrals write(f, "$ndihedrals dihedrals") end
        if hasimpropers write(f, "$nimpropers impropers") end
        write(f, "$natom_types atom types")
        if hasbonds write(f, "$nbond_types bond types") end
        if hasangles write(f, "$nangle_types angle types") end
        if hasdihedrals write(f, "$ndihedral_types dihedral types") end
        if hasimpropers write(f, "$nimproper_types improper types") end
        box = snapshot["box"]
        for (i, d) in enumerate(['x', 'y', 'z'])
            write(f, "$(box[i, 1]) $(box[i, 2]) $(d)lo $(d)hi\n")
        end
        if hastilt write(f, "$(snapshot["tilt"][1]) $(snapshot["tilt"][2]) $(snapshot["tilt"][3]) xy xz yz\n\n") end

        write(f, "Atoms # $atom_style\n\n")
        for i=1:natoms
            write(f, "$(snapshot["atomIDs"][i])")
            if hasbonds write(f, " $(snapshot["molecules"][i])") end
            write(f, " $(snapshot["atom_types"][i])")
            if hascharges write(f, " $(snapshot["charges"][i])") end
            write(f, " $(snapshot["coords"][i, 1]) $(snapshot["coords"][i, 2]) $(snapshot["coords"][i, 3])")
            if hasimages write(f, " $(snapshot["images"][i, 1]) $(snapshot["images"][i, 2]) $(snapshot["images"][i, 3])") end
            write(f, "\n")
        end

        if hasvelocities
            write(f, "\nVelocities\n\n")
            for (id, vel) in zip(snapshot["atomIDs"], eachrow(snapshot["velocities"]))
                write(f, "$id $(vel[1]) $(vel[2]) $(vel[3])\n")
            end
        end

        if hasbonds
            write(f, "\nBonds\n\n")
            for (i, (bt, bond)) in enumerate(zip(snapshot["bond_types"], eachrow(snapshot["bonds"])))
                write(f, "$i $bt $(bond[1]) $(bond[2])\n")
            end
        end

        if hasangles
            write(f, "\nAngles\n\n")
            for (i, (at, angle)) in enumerate(zip(snapshot["angle_types"], eachrow(snapshot["angles"])))
                write(f, "$i $at $(angle[1]) $(angle[2]) $(angle[3])\n")
            end
        end

        if hasdihedrals
            write(f, "\nDihedrals\n\n")
            for (i, (dt, dihedral)) in enumerate(zip(snapshot["dihedral_types"], eachrow(snapshot["dihedrals"])))
                write(f, "$i $dt $(dihedral[1]) $(dihedral[2]) $(dihedral[3]) $(dihedral[4])\n")
            end
        end

        if hasimpropers
            write(f, "\nImpropers\n\n")
            for (i, (bt, improper)) in enumerate(zip(snapshot["improper_types"], eachrow(snapshot["impropers"])))
                write(f, "$i $bt $(improper[1]) $(improper[2]) $(improper[3]) $(improper[4])\n")
            end
        end
    end
    return nothing
end
