"""
    Takes raw mathematica output and formats it to something that is way nicer for saving.
    string: Raw mathematica string
    args: List of how the χs are called, e.g. ["u1","d1","e1","u2"]
"""
function clean_string(string, args)
    string = replace(string, "\n" => "")
    string = filter(x -> !isspace(x), string)
    for arg in args
        string = replace_subscript!(string, arg)
    end
    return string
end

"""
    Helper for clean_string. Manipulates the "u1" or similar subscripts.
"""
function replace_subscript!(string, sub)
    string = replace(string, "Subscript[\\[Chi],"*sub*"]" => sub)
    return string
end

using HDF5
"""
    Save a nicely formatted string for later reuse.
    Use EoverN as a_b instead of a/b to not mess up structure of file.
"""
function save_string(string, args, EoverN, type)
    filename = args2string(args)
    h5write("./data/modelzoo/"*filename*".h5", EoverN*"/"*type, string)
end

"""
    Helper for save_string. Takes ["u1","d1","e1","u2"] and returns "u1d1e1u2".
"""
function args2string(args)
    string = ""
    for arg in args
        string *= arg 
    end
    return string
end

model = ["u1","d1","e1","u2"]
model_str = args2string(model)

str = clean_string(str, ["u1","d1","e1","u2"])
save_string(str, ["u1","d1","e1","u2"], "-1_3", "Sols")

fullset = h5read("./data/modelzoo/"*model_str*".h5", "-1_3")
fullset["Subsets"]









str = raw""