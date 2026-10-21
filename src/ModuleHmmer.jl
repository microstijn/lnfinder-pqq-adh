module ModuleHmmer

using BioSequences
using FASTX
using DataFrames

export parse_hmmalign_results, calibrate_active_site, identify_model_coordinate

"""
    identify_model_coordinate(hmm_aligned_seq::String, target_pattern::String)

Helper to find the index of a specific pattern (e.g., "DVD") in a single reference 
sequence string that has already been mapped to HMM coordinates.
Useful for debugging specific sequences.
"""
function identify_model_coordinate(hmm_aligned_seq::String, target_pattern::String)
    # Strip insertions (lowercase) to match HMM coordinate system
    clean_seq = filter(c -> isuppercase(c) || c == '-', hmm_aligned_seq)
    
    # Find the range of the pattern
    rng = findfirst(target_pattern, clean_seq)
    if isnothing(rng)
        error("Pattern '$target_pattern' not found in the reference sequence.")
    end
    # Return the start index (1-based)
    return rng.start
end

"""
    calibrate_active_site(fasta_path::String)

Reads the first two sequences from an aligned FASTA (AFA) to find the model coordinate.
Assumes:
  Seq 1 = Lanthanide Reference (e.g. XoxF1) -> Expects D-x-D
  Seq 2 = Calcium Reference    (e.g. MxaF)  -> Expects D-x-[ATS]

Source for that asymmetry: Good et al. 2020, "Lanthanide-dependent alcohol dehydrogenases
require an essential aspartate residue for metal coordination and enzymatic function",
J. Biol. Chem. 295:8272-8284, doi:10.1074/jbc.RA120.013227. It reports crystal structures of
XoxF1 from M. extorquens AM1 in which the lanthanide is coordinated by Glu192 (bidentate),
Asn276 (monodentate), Asp318 (monodentate) and Asp320 (bidentate). Asp318 and Asp320 are the
D-x-D pair this function looks for, and Good et al. report that nearly all XoxF and ExaF
sequences carry that ADDITIONAL aspartate at the position where the calcium enzymes carry
Ala, Ser or Thr (the D-x-[ATS] above). The pair is what makes the aligned column a Ca/Ln
switch rather than a conserved site, which is the only property this function relies on.

Returns the 1-based Model Coordinate (index ignoring insertions).
"""
function calibrate_active_site(fasta_path::String)
    
    # Read the first two sequences
    seqs = String[]
    reader = FASTA.Reader(open(fasta_path))
    
    for (i, record) in enumerate(reader)
        push!(seqs, FASTA.sequence(String, record))
        if i >= 2; break; end
    end
    close(reader)

    if length(seqs) < 2
        error("File contains fewer than 2 sequences. Cannot calibrate.")
    end

    # 2. Convert to Model Coordinates (Strip lowercase insertions)
    # This aligns them strictly to the HMM node numbers (1, 2, 3...)
    ref_ln = filter(c -> isuppercase(c) || c == '-', seqs[1])
    ref_ca = filter(c -> isuppercase(c) || c == '-', seqs[2])

    # Sanity Check: Model coordinates must be identical length for all sequences
    if length(ref_ln) != length(ref_ca)
        error("HMM Alignment error: Model coordinate strings are different lengths!")
    end

    println("Scanning for Active Site Collision...")
    println("Ref 1 (Ln): Looking for D-x-D")
    println("Ref 2 (Ca): Looking for D-x-[ATS]")
    println("-"^40)

    # 3. Scan for the unique overlap
    # Iterate through the MODEL string
    for i in 1:(length(ref_ln) - 2)
        # Check Reference 1 for Lanthanide Pattern (D-x-D)
        r1_anchor = ref_ln[i]
        r1_switch = ref_ln[i+2]
        is_ln   = (r1_anchor == 'D' && r1_switch == 'D')

        # Check Reference 2 for Calcium Pattern (D-x-A/T/S)
        r2_anchor = ref_ca[i]
        r2_switch = ref_ca[i+2]
        is_ca   = (r2_anchor == 'D' && (r2_switch in ['A', 'T', 'S']))

        # Check for COLLISION (The confirmation)
        if is_ln && is_ca
            println(">>> MATCH FOUND AT MODEL INDEX: $i")
            println("    Ref 1: $r1_anchor-$(ref_ln[i+1])-$r1_switch")
            println("    Ref 2: $r2_anchor-$(ref_ca[i+1])-$r2_switch")
            return i
        end
    end

    error("Could not find a column where Ref 1 is D-x-D and Ref 2 is D-x-[ATS]. Check your references.")
end

"""
    parse_hmmalign_results(fasta_path::AbstractString, active_site_idx::Int)

Parses the output of `hmmalign --trim --outformat afa`.
Captures the EXACT 3-residue motif for visualization.
"""
function parse_hmmalign_results(fasta_path::AbstractString, active_site_idx::Int)
    ids = String[]
    accessions = String[] # the merge key against the taxonomy and .uc tables
    motifs = Symbol[] 
    residues = String[] 

    FASTA.Reader(open(fasta_path)) do reader
        for record in reader
            raw_header = FASTA.identifier(record)
            
            acc = ""
            if occursin('|', raw_header)
                parts = split(raw_header, '|')
                acc = length(parts) >= 2 ? String(parts[2]) : String(parts[1])
            else
                m = match(r"^([A-Z0-9]+)", raw_header)
                acc = isnothing(m) ? raw_header : String(m.captures[1])
            end
            acc = split(acc, ' ')[1] 

            raw_seq = FASTA.sequence(String, record)
            model_seq = filter(c -> isuppercase(c) || c == '-', raw_seq)
            
            if active_site_idx + 2 > length(model_seq)
                push!(ids, raw_header)
                push!(accessions, acc)
                push!(motifs, :Truncated)
                push!(residues, "---")
                continue
            end

            actual_motif = model_seq[active_site_idx : active_site_idx + 2]
            anchor = actual_motif[1] 
            switch = actual_motif[3] 
            
            if anchor == '-' || switch == '-' || actual_motif[2] == '-' 
                motif_type = :Gap
            elseif anchor == 'D'
                if switch == 'D'
                    motif_type = :Lanthanide
                elseif switch in ['A', 'T', 'S']
                    motif_type = :Calcium
                else
                    motif_type = :Unknown 
                end
            else
                motif_type = :Unknown
            end

            push!(ids, raw_header)
            push!(accessions, acc)
            push!(motifs, motif_type)
            push!(residues, String(actual_motif)) 
        end
    end

    return DataFrame(original_header=ids, accession=accessions, classification=motifs, active_site=residues)
end

end # module