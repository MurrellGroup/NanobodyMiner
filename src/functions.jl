#Copyright: Ben Murrell
#All rights reserved.
#A license clarifying use of this code has not yet been provided

#top level function 

#This should be moved into seqUMAP
function seqPCA(seqs::Array{String,1}; 
    k = 5,
    lookup_dic = NT_DICT,
    pca_maxoutdim = 5)
    vecs = []
    missing_chars = Set{Char}()
    for seq in seqs
        vec, missing_chars = SeqUMAP.kmer_embed(seq, k, kmer_count!; lookup_dic = lookup_dic, missing_chars = missing_chars);
        push!(vecs, vec)
    end
    if length(missing_chars) > 0
        @warn "Ignored the following characters missing from lookup: $(unique(missing_chars)). Check your sequence type!"
    end
    X = hcat(vecs...);
    X = convert(Array{Float64,2}, X)

    M = fit(SeqUMAP.PCA, X; maxoutdim = pca_maxoutdim)
    PCAembedding = SeqUMAP.transform(M, X);        
    return PCAembedding
end

function collapse(seqs; prefix = "")
    count_dic = countmap(seqs);
    keep_arr = sort([(count_dic[k],k) for k in keys(count_dic) if count_dic[k]>0],rev = true)
    return [prefix*"$(i)_$(ka[1])" for (i,ka) in enumerate(keep_arr)],[ka[2] for (i,ka) in enumerate(keep_arr)]
end

robust_translate(nuc_seq) = translate_to_aa(nuc_seq[1:3*Int(floor(end/3))])
size_from_name(n) = parse(Int64,split(n,"_")[end])
dataset_from_name(n) = split(n,"_")[1]
#barcode_from_prefix(p; d = barcode2int) = get(d,p[[4,19,21,24]],0)
#Hardcoded positions of bases to use as barcodes.
#If error, then switch back to this one:
#barcode_from_prefix(p; d = barcode2int) = get(d,p[[4,19,21]],0)
barcode_from_prefix(p; d = barcode2int, inds = [4,19,21]) = get(d,p[inds],0)
hamming_d(s1,s2) = sum(collect(s1) .!= collect(s2))

#takes in an alignment, a gap character, and returns:
#a dictionary that maps the degapped alignment[1] coordinates to the degapped alignment[2] coordinates
#Coordinates are returned as tuples, because:
#If a position in alignment[1] maps to a gap in alignment[2],
#then the retuned tuple coords are the query coords either side of the gap
#If a position in alignment[1] maps to a non-gap in alignment[2], the the tuple values are equal.
#. 12 3456789
#1:GT-ACTWAGT
#2:GTWACT-A-T
#  123456 7 8
#Edge case: the map(6) and map(8) map to gaps, but sometimes we want map(6)=7 and sometimes we want map(8)=7
#Eg. when we want to pull out a CDR region.
#So we need to know whether we're pulling out a "start" or an "end" coord.
#So when eg. map(x) maps to a gap, we return the inds either side of the gap, so the user can decide.
#In the above example, map(6) = (6,7), map(7) = (7,7), and map(8) = (7,8).
function region_map(ali; gap = '-')
    map_dict = Dict{Int64, Tuple{Int64, Int64}}()
    ref_coord = 1
    query_coord = 1
    for i in 1:length(ali[1])
        refgap = ali[1][i] == gap
        querygap = ali[2][i] == gap
        if refgap
            if !querygap
                query_coord += 1
            end
            continue
        end
        if !refgap && querygap
            map_dict[ref_coord] = (query_coord-1, query_coord)
            ref_coord += 1
            continue
        end
        if !refgap && !querygap
            map_dict[ref_coord] = (query_coord, query_coord)
            ref_coord += 1
            query_coord += 1
        end
    end
    return map_dict
end

#  12 345678910
#R:GT-ACTWAGGT
#Q:GTWACT-A--T
#  123456 7  8
#sort(region_map(["GT-ACTWAGGT",
#                 "GTWACT-A--T"]))

#These CDRs need to be specific to the consensus sequence from the gaptrimmed antibody alignment
#Returns the coordinates in the original sequence of the regions after profile alignment to the nanobody profile
#I should pull these coords from the top line of the nanobody MSA?
function CDRcut(s; coords = [(19,27),(39,51),(89,107),(1,118)], refprof = nanobodyprofile)
    ali = profile_affine_align(nanobodyprofile,seqs2profile([s]), profile_cost, gap_extend = -0.9)
    rm = region_map(ali, gap = [('!', 1.0)])
    region_coords = [(rm[c[1]][2],rm[c[2]][1]) for c in coords]
    return region_coords
end

#=
#This is if you need to generate a new 
write_fasta("random1000.fasta",sample(df.trimmed_AA_seq,1000));
#Then align these, outside
gappy_AA_aligned = read_fasta("random1000_gappy_aligned.fasta")[2];
function simple_consensus(s_arr)
    [mode([s[i] for s in s_arr]) for i in 1:length(s_arr[1])];
end
col_keeps = [sum([s[i] == '-' for s in gappy_AA_aligned])/length(gappy_AA_aligned) < 0.5 for i in 1:length(gappy_AA_aligned[1])];
gaptrimmed = [join(collect(g)[col_keeps]) for g in gappy_AA_aligned];
write_fasta("random1000_gaptrimmed.fasta",gaptrimmed);

#Then run this to get the consensus, and use that to get CDR coords
s_arr = gaptrimmed
join([mode([s[i] for s in s_arr if s[i] != '-']) for i in 1:length(s_arr[1])])
=#

function extract_AA_regions(AA,regions)
    return [AA[r[1]:r[2]] for r in regions]
end

#a tuple of start and end AA coords
function aa2nuc(coord)
    return 1+3*(coord[1]-1):3+3*(coord[2]-1)
end

function extract_nuc_regions(nuc,AA_region_coords)
    return [nuc[aa2nuc(r)] for r in AA_region_coords]
end



#Functions for selecting candidates that are different from each other
function robust_hamming_prop(s1,s2)
    if length(s1) != length(s2)
        return 1.0
    else
        return Distances.hamming(s1,s2)/length(s1)
    end
end

function edit_dist_ratio(s1,s2)
    kmer_seeded_edit_dist(s1,s2) ./ ((length(s1)+length(s2))/2)
end

function prioritize(scores, number; 
        lineages = nothing, 
        previously_selected = Int64[], 
        AAs = nothing, 
        AAthresh = 0.15,
        CDR3s = nothing,
        CDR3thresh = 0.4)
    selected_set = copy(previously_selected)
    score_rank = sortperm(scores, rev = true)
    i = 1
    for i in 1:10000
        if length(selected_set) >= number
            break
        end
        ind = score_rank[i]
        
        lineage_flag = true
        if !isnothing(lineages)
            if lineages[ind] in lineages[selected_set]
                lineage_flag = false
            end
        end
        
        CDR3_flag = true
        if lineage_flag #This just skips a bunch of alignments when its already rules out
            if (!isnothing(CDR3s)) && (length(selected_set) >= 1)
                if minimum([robust_hamming_prop(CDR3s[ind],s) for s in CDR3s[selected_set]]) < CDR3thresh
                    CDR3_flag = false
                end
            end
        end
        
        ham_flag = true
        if lineage_flag && CDR3_flag #This just skips a bunch of alignments when its already rules out
            if (!isnothing(AAs)) && (length(selected_set) >= 1)
                if minimum([edit_dist_ratio(AAs[ind],s) for s in AAs[selected_set]]) < AAthresh
                    ham_flag = false
                end
            end
        end
        if lineage_flag && CDR3_flag && ham_flag
            push!(selected_set,ind)
        end
    end
    return selected_set
end

#For matching one element in a list, picking the largest match
function get_biggest_ind(query,list,sizes)
    inds = findall(list .== query)
    if length(inds) > 0
        return inds[argmax(sizes[inds])]
    else
        return nothing
    end
end

#Looks for a perfect or partial seq match. Returns nothing, or the match index and distance if a match
function fast_search(quer,AAlist,sizes; cutoff = 10)
    found_ind = get_biggest_ind(quer,AAlist,sizes)
    if !isnothing(found_ind)
        return found_ind, 0
    else
        dists = [kmer_seeded_edit_dist(s,quer) for s in AAlist]
        mini = argmin(dists)
        if dists[mini] < cutoff
            inds = findall(dists[mini] .== dists)
            if length(inds) > 0
                return inds[argmax(sizes[inds])],dists[mini]
            else
                return nothing
            end
        end
    end
end

#This stuff handles the implicit "barcoding" with the degenerate bases in the cloning primers
function grow(arr::Array{String},elems::Array{String})
    new_arr = String[]
    for a in arr
        for elem in elems
            push!(new_arr,a*elem)
        end
    end
    return new_arr
end

#This needs to be made more robust
#returns dataset,count,barcode_num,prefix,suffix, hinge,functional,seq,AAseq
function encode_variant(name,seq)
    
    prefix_match = "CAGNTGCAGCTCGTGGAGNCNGGN"
    
    long_hinge_suffix  = "CCCAAGACACCAAAACCACAACCGGCGCGCCAGGCC"
    short_hinge_suffix = "GCGCACCACAGCGAAGACCCCTCGGCGCGCCAGGCC"
    suffix_length = length(short_hinge_suffix)
    
    count = size_from_name(name)
    dataset = dataset_from_name(name)
    prefix = seq[1:24] #This is the last codon with an ambig.
    barcode_num = barcode_from_prefix(prefix)
    AA_seq = robust_translate(seq[25:end])
    functional = true
    if AA_seq[end-3:end] != "ARQA" || occursin("*",AA_seq) || hamming_d(prefix,prefix_match) > 6
        functional = false
    end
    
    suffix = seq[end+1-suffix_length:end]
    
    dists = [hamming_d(suffix,long_hinge_suffix),hamming_d(suffix,short_hinge_suffix)]
    hinge = "bad"
    if minimum(dists) < 2
        hinge = ["long","short"][argmin(dists)]
    end
    
    if hinge == "bad"
        outseq = seq[26:end]
        elseif hinge == "short"
            outseq = seq[25:end-suffix_length]
        else
            outseq = seq[25:end-(suffix_length+3)]
    end
    trimmed_AA_seq = robust_translate(outseq)
    
    return dataset,count, barcode_num,prefix,hinge,functional,trimmed_AA_seq,outseq
end

#=
#This is the consensus of the current panel nanobody.
>ConsensusWithPrefix
QVQLVESGGGLVQPGGSLRLSCAASGFTFDDYAIGWFRQAPGKEREGVSCISSSDGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAERGGYYSDCITGYGYDYWGQGTQVTVSS

AbNum Results, Chothia
H1 Q
H2 V
H3 Q
H4 L
H5 V
H6 E
H7 S
H8 G
H9 G
H10 G
H11 L
H12 V
H13 Q
H14 P
H15 G
H16 G
H17 S
H18 L
H19 R
H20 L
H21 S
H22 C
H23 A
H24 A
H25 S
H26 G
H27 F
H28 T
H29 F
H30 D
H31 D
H32 Y
H33 A
H34 I
H35 G
H36 W
H37 F
H38 R
H39 Q
H40 A
H41 P
H42 G
H43 K
H44 E
H45 R
H46 E
H47 G
H48 V
H49 S
H50 C
H51 I
H52 S
H52A S
H53 S
H54 D
H55 G
H56 S
H57 T
H58 Y
H59 Y
H60 A
H61 D
H62 S
H63 V
H64 K
H65 G
H66 R
H67 F
H68 T
H69 I
H70 S
H71 R
H72 D
H73 N
H74 A
H75 K
H76 N
H77 T
H78 V
H79 Y
H80 L
H81 Q
H82 M
H82A N
H82B S
H82C L
H83 K
H84 P
H85 E
H86 D
H87 T
H88 A
H89 V
H90 Y
H91 Y
H92 C
H93 A
H94 A
H95 E
H96 R
H97 G
H98 G
H99 Y
H100 Y
H100A S
H100B D
H100C C
H100D I
H100E T
H100F G
H100G Y
H100H G
H100I Y
H101 D
H102 Y
H103 W
H104 G
H105 Q
H106 G
H107 T
H108 Q
H109 V
H110 T
H111 V
H112 S
H113 S
=#