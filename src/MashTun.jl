module MashTun

export minhash,
       jaccarddist,
       mashdist

using Bio.Seq
using DataStructures


function rchash(kmer::Kmer)
    k = minimum((kmer, reverse_complement(kmer)))
    return hash(k)
end


function kmerminhash(seq::BioSequence, kmerset, kmerhashes::Vector{UInt64}, k::Int, s::Int)
    typeof(kmerset) <: Set || typeof(kmerset) <: SortedSet ? true : error("Kmerset must be a `Set` or `SortedSet`")

    for kmer in each(DNAKmer{k}, seq)
        if length(kmerhashes) == 0
            if length(kmerset) < s
                push!(kmerset, hash(kmer[2]))
            elseif length(kmerset) == s
                kmerset = SortedSet(kmerset)
                for i in kmerset
                    push!(kmerhashes, pop!(kmerset))
                end
            end
        end

        if length(kmerhashes) == s
            h = hash(kmer[2])
            if h < kmerhashes[end]
                i = searchsortedlast(kmerhashes, h)
                if i == 0 && i != kmerhashes[1]
                    pop!(kmerhashes)
                    unshift!(kmerhashes, h)
                elseif h != kmerhashes[i]
                    pop!(kmerhashes)
                    insert!(kmerhashes, i+1, h)
                end
            end
        end

    end
    return (kmerset, kmerhashes)
end


function minhash(seq::BioSequence, k::Int, s::Int)
    kmerset = Set{UInt64}()
    kmerhashes = Vector{UInt64}()
    kmerset, kmerhashes = kmerminhash(seq, kmerset, kmerhashes, k, s)
    return kmerhashes
end

function minhash{T<:BioSequence}(seqs::Vector{T}, k::Int, s::Int)
    kmerset = Set{UInt64}()
    kmerhashes = Vector{UInt64}()

    for seq in seqs
        kmerset, kmerhashes = kmerminhash(seq, kmerset, kmerhashes, k, s)
    end
    return kmerhashes
end


function minhash{T<:BioSequence}(seqs::FASTAReader{T}, k::Int, s::Int)
    kmerset = Set{UInt64}()
    kmerhashes = Vector{UInt64}()
    for seq in seqs
        kmerset, kmerhashes = kmerminhash(seq.seq, kmerset, kmerhashes, k, s)
    end
    return kmerhashes
end


function jaccarddist(sketch1::Vector{UInt64}, sketch2::Vector{UInt64}, s::Int)
    i = 0
    matches = 0
    sk1 = deepcopy(sketch1)
    sk2 = deepcopy(sketch2)

    n1 = shift!(sk1)
    n2 = shift!(sk2)

    while i < s && length(sk1) != 0 && length(sk2) != 0
        if n1 == n2
            matches += 1
            i += 1
            n1 = shift!(sk1)
            n2 = shift!(sk2)
        elseif n1 < n2
            while n1 < n2
                i += 1
                n1 = shift!(sk1)
            end
        elseif n2 < n1
            while n2 < n1
                i += 1
                n2 = shift!(sk2)
            end
        end
    end
    println(matches)
    return matches / i
end


function mashdist(k::Int, j::Float64)
    return 1/k * log(2j / (1+j))
end


end # module MashTun
