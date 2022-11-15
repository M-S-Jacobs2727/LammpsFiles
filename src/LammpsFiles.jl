module LammpsFiles
export readData, readDump, wrap, unwrap, writeData
# TODO: unwrap, write_data, read_data
using Logging

include("dump.jl")

include("data.jl")

include("wrap.jl")

end
