using Logging

struct snapshot
    atom_ids::Vector{<:Integer}
    idtoindex::Vector{<:Integer}
    molecules::Vector{<:Integer}
    atom_types::Vector{<:Integer}
    charges::Vector{<:Real}
    coords::Matrix{<:Real}
    images::Matrix{<:Integer}
    velocities::Matrix{<:Real}
    the_rest::Matrix{<:Real}
    bonds::Matrix{<:Integer}
    angles::Matrix{<:Integer}
    dihedrals::Matrix{<:Integer}
    impropers::Matrix{<:Integer}
    bond_types::Vector{<:Integer}
    angle_types::Vector{<:Integer}
    dihedral_types::Vector{<:Integer}
    improper_types::Vector{<:Integer}
end

function removeComment(line)
    strip(split(line, "#")[1])
end

function readDataSection(f, data)
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
function readData(source, atom_style="full")
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

    open(source, "r") do f
        readline(f)  # skip first line
        line = readline(f)  # next line blank

        # Header values
        natoms = natom_types = nbonds = nbond_types = nangles = nangle_types = \ 
            ndihedrals = ndihedral_types = nimpropers = nimproper_types = 0
        box = zeros(3, 2)
        tilt = zeros(3)

        image_flag = false

        # Parse header, read until non-blank line without header keyword
        while true
            line = readline(f)
            line = removeComment(line)
            if line == ""
                continue
            end
            
            if occursin("atoms", line)
                natoms = parse(Int, split(line, " atoms")[1])
            elseif occursin("bonds", line)
                nbonds = parse(Int, split(line, " bonds")[1])
            elseif occursin("angles", line)
                nangles = parse(Int, split(line, " angles")[1])
            elseif occursin("dihedrals", line)
                ndihedrals = parse(Int, split(line, " dihedrals")[1])
            elseif occursin("impropers", line)
                nimpropers = parse(Int, split(line, " impropers")[1])
            elseif occursin("atom types", line)
                natom_types = parse(Int, split(line, " atom types")[1])
            elseif occursin("bond types", line)
                nbond_types = parse(Int, split(line, " bond types")[1])
            elseif occursin("angle types", line)
                nangle_types = parse(Int, split(line, " angle types")[1])
            elseif occursin("dihedral types", line)
                ndihedral_types = parse(Int, split(line, " dihedral types")[1])
            elseif occursin("improper types", line)
                nimproper_types = parse(Int, split(line, " improper types")[1])
            elseif occursin("xlo xhi", line)
                values = split(line)
                box[1, : ] = [parse(Float32, values[1]), parse(Float32, values[2])]
            elseif occursin("ylo yhi", line)
                values = split(line)
                box[1, : ] = [parse(Float32, values[1]), parse(Float32, values[2])]
            elseif occursin("zlo zhi", line)
                values = split(line)
                box[1, : ] = [parse(Float32, values[1]), parse(Float32, values[2])]
            elseif occursin("xy xz yz", line)
                values = split(line)
                tilt = [parse(Float32, values[1]), parse(Float32, values[2]), parse(Float32, values[3])]
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

        # Body values
        atom_ids = zeros(Int, natoms)
        molecules = zeros(Int, natoms)
        atom_types = zeros(Int, natoms)
        charges = zeros(Float32, natoms)
        coords = zeros(Float32, natoms, 3)
        images = zeros(Int, natoms, 3)
        velocities = zeros(Float32, natoms, 3)
        the_rest = zeros(Float32, natoms, 2)

        # Molecule body values
        bonds = zeros(Int, nbonds, 2)
        angles = zeros(Int, nangles, 3)
        dihedrals = zeros(Int, ndihedrals, 4)
        impropers = zeros(Int, nimpropers, 4)
        bond_types = zeros(Int, nbond_types)
        angle_types = zeros(Int, nangle_types)
        dihedral_types = zeros(Int, ndihedral_types)
        improper_types = zeros(Int, nimproper_types)

        # Body sections supported: Atoms, Velocities, Masses, Bonds, Angles, Dihedrals, Impropers,
        #       Pair Coeffs, PairIJ Coeffs, Bond Coeffs, Angle Coeffs, Dihedral Coeffs, Improper Coeffs
        # Body sections
        indices = zeros(Int, natoms)  # for sorting
        while true
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
                first_row_data = parse.(Float32, split(removeComment(line)))
                
                num_cols = length(first_row_data)
                if num_cols != base_num_cols && num_cols != base_num_cols + 3
                    throw(DimensionMismatch("invalid number of columns in Atoms section of data file: \
                                         $(num_cols). Expected $(base_num_cols) or $(base_num_cols+3)"))
                end
                if num_cols == base_num_cols + 3
                    image_flag = true
                end

                tmp_data = zeros(num_cols, natoms)
                tmp_data[ : , 1] = first_row_data
                tmp_data[ : , 2:end] = readDataSection(f, tmp_data[ : , 2:end])

                indices = sortperm(tmp_data[1, : ])
                tmp_data = tmp_data[ : , indices]
                
                atom_ids = Integer.(tmp_data[1, : ])
                if atom_style == "atomic"
                    atom_types = Integer.(tmp_data[2, : ])
                    coords = tmp_data[3:5, : ]'
                elseif atom_style == "charge"
                    atom_types = Integer.(tmp_data[2, : ])
                    charges = tmp_data[3, : ]
                    coords = tmp_data[4:6, : ]'
                elseif (
                    atom_style == "molecular"
                    || atom_style == "bond"
                    || atom_style == "angle"
                )
                    molecules = Integer.(tmp_data[2, : ])
                    atom_types = Integer.(tmp_data[3, : ])
                    coords = tmp_data[4:6, : ]'
                elseif atom_style == "full"
                    molecules = Integer.(tmp_data[2, : ])
                    atom_types = Integer.(tmp_data[3, : ])
                    charges = tmp_data[4, : ]
                    coords = tmp_data[5:7, : ]'
                elseif (
                    atom_style == "line
                    mesont
                    bpm/sphere
                    tri"
                )
                    molecules = Integer.(tmp_data[2, : ])
                    atom_types = Integer.(tmp_data[3, : ])
                    the_rest = image_flag ? tmp_data[4:end-6, : ]' : tmp_data[4:end-3, : ]'
                    coords = image_flag ? tmp_data[end-5:end-3, : ]' : tmp_data[end-2:end, : ]'
                else
                    atom_types = Integer.(tmp_data[2, : ])
                    the_rest = image_flag ? tmp_data[4:end-6, : ]' : tmp_data[4:end-3, : ]'
                    coords = image_flag ? tmp_data[end-5:end-3, : ]' : tmp_data[end-2:end, : ]'
                end
                if image_flag
                    images = Integer.(tmp_data[end-2:end, : ]')
                end
                idtoindex = zeros(maximum(atom_ids))
                for (i, a) in enumerate(atom_ids)
                    idtoindex[a] = i
                end

            elseif occursin("Velocities", line)
                readline(f)  # Blank
                data = zeros(4, natoms)
                data = readDataSection(f, data)
                # Try applying indices from Atoms section to sort velocities by ID
                # Check if sorted. If not, sort normally. Will usually save time.
                issorted(data[1, indices]) || (indices = sortperm(data[1, : ]))
                velocities = data[2:4, indices]'
            elseif occursin("Masses", line)  # Not yet implemented
                readline(f)  # Blank
                for i in 1:natom_types; readline(f); end
            elseif nbonds > 0 && occursin("Bonds", line)
                readline(f)  # Blank
                data = zeros(Integer, 4, nbonds)
                data = readDataSection(f, data)
                bond_types = data[2, : ]
                bonds = data[3:4, : ]'
            elseif nangles > 0 && occursin("Angles", line)
                readline(f)  # Blank
                data = zeros(Integer, 5, nangles)
                data = readDataSection(f, data)
                angle_types = data[2, : ]
                angles = data[3:5, : ]'
            elseif ndihedrals > 0 && occursin("Dihedrals", line)
                readline(f)  # Blank
                data = zeros(Integer, 6, ndihedrals)
                data = readDataSection(f, data)
                dihedral_types = data[2, : ]
                dihedrals = data[3:6, : ]'
            elseif nimpropers > 0 && occursin("Impropers", line)
                readline(f)  # Blank
                data = zeros(Integer, 7, nimpropers)
                data = readDataSection(f, data)
                improper_types = data[2, : ]
                impropers = data[3:6, : ]'
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
            if eof(f)
                break
            end
        end
    end
    snapshot(
        atom_ids,
        idtoindex,
        molecules,
        atom_types,
        charges,
        coords,
        images,
        velocities,
        the_rest,
        bonds,
        angles,
        dihedrals,
        impropers,
        bond_types,
        angle_types,
        dihedral_types,
        improper_types
    )
end

function writeData(source, snapshot, atom_style)
    return nothing
end
