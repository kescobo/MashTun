module MashTun

export minhash,
       mashdistance,
       MASHSketch

using Bio.Seq

include("src/minhash.jl")
include("src/mash.jl")

end # module MashTun
