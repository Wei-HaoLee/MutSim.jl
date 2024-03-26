using FASTX

function read_fa(path::String)
    rdr = FASTAReader(open(path, "r"));
    records = collect(rdr); close(rdr);

    return records
end

function read_fq(path::String)
    rdr = FASTQReader(open(path, "r"))
    records = collect(rdr); close(rdr);

    return records
end

function read_seq(seq::String)
end

function get_seq(record::Vector{FASTX.FASTA.Record})
    return map(rec -> sequence(String, rec), record)
end

function get_description(record::Vector{FASTX.FASTA.Record})
    return map(description, record)
end