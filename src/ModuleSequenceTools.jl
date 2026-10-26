
module SequenceTools

using BioSequences
using DataFrames
using FASTX

export count_residues_per_position, export_tree_fasta

"""
    count_residues_per_position(motif_column::Vector{String})
    
Returns a DataFrame where rows are Amino Acids and columns are 
the probability/count at each of the 3 positions in the D-x-D motif.
"""
function count_residues_per_position(motif_column)
    # 1. Filter out any remaining gaps or truncated entries
    valid_motifs = [m for m in motif_column if length(m) == 3 && !occursin("-", m)]
    total = length(valid_motifs)
    
    # 2. Get the standard amino acid alphabet
    amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    
    # 3. Initialize counts
    counts = Dict(aa => [0, 0, 0] for aa in amino_acids)
    
    # 4. Populate counts
    for motif in valid_motifs
        for pos in 1:3
            aa = string(motif[pos])
            if haskey(counts, aa)
                counts[aa][pos] += 1
            end
        end
    end
    
    # 5. Convert to DataFrame
    df_counts = DataFrame(AminoAcid = amino_acids)
    df_counts.Pos1 = [counts[aa][1] for aa in amino_acids]
    df_counts.Pos2 = [counts[aa][2] for aa in amino_acids]
    df_counts.Pos3 = [counts[aa][3] for aa in amino_acids]
    
    # 6. Add percentages alongside the raw counts
    df_counts.Pct1 = round.((df_counts.Pos1 ./ total) .* 100, digits=2)
    df_counts.Pct2 = round.((df_counts.Pos2 ./ total) .* 100, digits=2)
    df_counts.Pct3 = round.((df_counts.Pos3 ./ total) .* 100, digits=2)
    
    return sort(df_counts, :Pos2, rev=true)
end


"""
    export_tree_fasta(df_dedup::DataFrame, input_fasta::String, output_fasta::String)
Streams sequences from disk and writes a FASTA whose headers carry metadata
(`id|kingdom|habitat|pos2`), with explicit newlines between records.

!!! warning "Legacy header format; kept for src/scriptLn.jl only"
    The headers written here are `id|kingdom|habitat|pos2`. A pipe-separated label is
    iTOL LCA syntax — it names an internal node by two of its descendants — so such a
    label is not a leaf id, and iTOL datasets keyed by accession will not bind to a tree
    built from this FASTA. The habitat field is also the caller's `report_group` column
    as it stands — whichever habitat aggregation that DataFrame carries — baked
    irreversibly into the leaf name.

    For a tree that has to bind to iTOL datasets, write plain `>ACCESSION` headers and
    carry every annotation in separate iTOL dataset files instead.

    These headers, not bare accessions, are what a tree built from this FASTA carries on its
    leaves. A Newick reader has to keep `_` verbatim for the spaces this function substitutes
    out of the habitat field to survive parse → prune → root → write; it is the pipe, not the
    underscore, that makes the label something other than a leaf id.
"""
function export_tree_fasta(df_dedup::DataFrame, input_fasta::String, output_fasta::String)
    # IDs to keep (cluster centroids)
    target_ids = Set(df_dedup.id)
    
    # Map IDs to header metadata
    metadata = Dict(row.id => (
        kingdom = coalesce(row.kingdom, "Unknown"),
        habitat = replace(coalesce(row.report_group, "Other"), " " => "_"),
        pos2 = row.pos2_aa
    ) for row in eachrow(df_dedup))

    println("🧵 Streaming sequences from $input_fasta...")
    
    count = 0
    # Read the original FASTA with FASTX
    FASTA.Reader(open(input_fasta)) do reader
        # Open the output file for writing
        open(output_fasta, "w") do writer
            for record in reader
                id = FASTA.identifier(record)
                
                if id in target_ids
                    meta = metadata[id]
                    # New header: ID|Kingdom|Habitat|Pos2
                    new_header = "$(id)|$(meta.kingdom)|$(meta.habitat)|$(meta.pos2)"
                    
                    # Sequence of the current record
                    seq = FASTA.sequence(record)
                    
                    # Write the record manually to control line endings
                    write(writer, ">", new_header, "\n")
                    write(writer, seq, "\n")
                    
                    count += 1
                end
            end
        end
    end
    println("✅ Export complete: wrote $count centroids to $output_fasta")
end

end