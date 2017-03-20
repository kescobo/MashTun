module MashTun

export minhash,
       mashdistance,
       MinHashSketch

using Bio.Seq

include("minhash.jl")
include("mash.jl")

end # module MashTun
