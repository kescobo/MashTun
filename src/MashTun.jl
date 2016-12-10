module MashTun

export minhash,
       mashdist

using Bio.Seq
using Logging
using DataStructures

function rchash(kmer::Kmer)
    k = minimum((kmer, reverse_complement(kmer)))
    return hash(k)
end

function gethashheap(sequence::BioSequence, k::Int)
    hashheap = binary_minheap(UInt64)
    for kmer in each(DNAKmer{k}, sequence)
        h = rchash(kmer[2])
        push!(hashheap, h)
    end
    return hashheap
end

function gethashheap{T<:BioSequence}(sequences::Vector{T}, k::Int)
    hashheap = binary_minheap(UInt64)

    for sequence in sequences
        for kmer in each(DNAKmer{k}, sequence)
            h = rchash(kmer[2])
            push!(hashheap, h)
        end
    end
    return hashheap
end

function gethashheap{T<:BioSequence}(sequences::FASTAReader{T}, k::Int)
    hashheap = binary_minheap(UInt64)

    for sequence in sequences
        for kmer in each(DNAKmer{k}, sequence.seq)
            h = rchash(kmer[2])
            push!(hashheap, h)
        end
    end
    return hashheap
end


function minhash(hashheap::BinaryHeap, k::Int, s::Int)
    last = nothing
    sketch = Deque{UInt64}()

    while length(sketch) < s
        h = pop!(hashheap)
        if h != last
            push!(sketch, h)
            last = h
        end
    end

    return sketch
end


function minhash(sequence::BioSequence, k::Int, s::Int)
    hashheap = gethashheap(sequence, k)
    minhash(hashheap, k, s)
end

function minhash{T<:BioSequence}(sequences::Vector{T}, k::Int, s::Int)
    hashheap = gethashheap(sequences, k)
    minhash(hashheap, k, s)
end

function minhash{T<:BioSequence}(sequences::FASTAReader{T}, k::Int, s::Int)
    hashheap = gethashheap(sequences, k)
    minhash(hashheap, k, s)
end



function jaccarddist(sketch1::Deque{UInt64}, sketch2::Deque{UInt64}, s::Int)
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
    return matches / i
end


function mashdist(k::Int, j::Float64)
    return 1/k * log(2j / (1+j))
end


end # module MashTun
