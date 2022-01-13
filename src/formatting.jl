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
    Use EoverN as "a_b" instead of a/b to not mess up structure of file.
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
        string *= arg*"_" 
    end
    string = string[1:end-1] # remove last "_"
    return string
end

"""
    Inverse of args2string
"""
function string2args(model_str)
    args = split(model_str, "_")
    return args
end


using DataFrames
"""
"""
function string2dataframe(model_file, EoverN)
    set = h5read("./data/modelzoo/"*model_file*".h5", EoverN)
    sols = set["Sols"]
    subsets = filleting_subsets(set, model_file)

    return subsets
end

function filleting_subsets(set, model_file)
    subsets = set["Subsets"][2:end-1]
    subsets = split(subsets, "},{")
    for i in 1:length(subsets)
        subsets[i] = replace(subsets[i], "{" => "")
        subsets[i] = replace(subsets[i], "}" => "")
    end
    subsets = split.(subsets, ",")

    model = string2args(model_file)
    # model_arr Index#1 -> nr of subset, Index#2 -> nr of condition, Index#3 -> nr of higgs, Higgs sortation same as in filename.
    #model_arr = zeros(length(subsets), length(subsets[1]), length(model))
    #= for subset in subsets
        for condition in subset
            for higgs in model
                ind = findfirst(higgs, condition)
                try
                    if condition[ind[1]-1] == "+"
                        println("1")
                    elseif condition[ind[1]-1] == "-"
                        println("-1")
                    elseif condition[ind[1]-1] !=
                        print
                    end
                catch
                    nothing
                end
            end
        end
     end=#
end


string2dataframe(model_str, "-1_3")





model = ["u1","d1","e1","u2"]
model_str = args2string(model)
model_str
a = findfirst("d2", model_str)
print(a)

str = clean_string(str, ["u1","d1","e1","u2"])
save_string(str, ["u1","d1","e1","u2"], "-1_3", "Sols")



fullset = h5read("./data/modelzoo/"*model_str*".h5", "-1_3")
fullset["Subsets"]









str = raw""