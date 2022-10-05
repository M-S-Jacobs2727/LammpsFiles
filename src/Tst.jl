module Tst
include("LammpsFiles.jl")
import .LammpsFiles

using Test

@testset "read_dump 2 atoms" begin
    frame = LammpsFiles.read_dump("../sample_dumps/2atoms.dump")
    @debug "$(frame.natoms)\n$(frame.timestep)\n$(frame.box)\n$(frame.properties)\n$(frame.atoms)"
    @test frame.natoms == 2
    @test frame.timestep == 0
    @test frame.box == Float32[-10 10; -10 10; -10 10]
    @test frame.properties == ["id", "mol", "type", "xu", "yu", "zu", "pxx", "pyy", "pzz"]
    @test frame.atoms == Float32[1.0 2.0; 0.0 0.0; 1.0 1.0; 0.56 0.57; 0.56 11.24; 0.56 0.88; 1.07 1.07; 1.07 1.07; 1.08 1.08]
end

@testset "wrap 2 atoms" begin
    frame = LammpsFiles.read_dump("../sample_dumps/2atoms.dump")
    @test isapprox(LammpsFiles.wrap(frame.atoms[4:6, : ], frame.box), Float32[0.56 0.57; 0.56 -8.76; 0.56 0.88])
end
end
