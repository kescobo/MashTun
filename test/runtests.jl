using Base.Test

using Bio.Seq
using MashTun


function random_seq(n::Integer, nts, probs)
    cumprobs = cumsum(probs)
    x = Array(Char, n)
    for i in 1:n
        x[i] = nts[searchsorted(cumprobs, rand()).start]
    end
    return convert(AbstractString, x)
end

function random_dna(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'T', 'N'], probs)
end

@testset "MashTun" begin
    @testset "MinHash" begin
        seq = DNASequence(random_dna(1000))
        h = minhash(seq, 10, 100)

        @test length(h) == 100
        @test h == minhash(seq, 10, 100)

        @test_throws BoundsError h[101]

    end
    @testset "MASH" begin

        a = minhash(dna"ATCGCCA-", 4, 3)
        b = minhash(dna"ATCGCCTA", 4, 3)
        @test_approx_eq_eps mashdistance(a, b) 0.2745 1e-3
        @test mashdistance(a, a) == 0 && mashdistance(b, b) == 0
        @test a.sketch == sort(a.sketch)
    end
end
