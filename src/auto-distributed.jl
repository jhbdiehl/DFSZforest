export start_of_shift, end_of_work

using Distributed

function start_of_shift(workers=length(Sys.cpu_info()))
    addprocs(workers)
    @everywhere println(pwd())
    println("Sent $workers workers to the factory.") 
end

function end_of_work()
    ids = workers()
    rmprocs(ids)
    num = length(ids)
    println("$num workers are going home to have a Feierabendbier.")
    if nprocs() != 1
        i = nprocs()
        println("You still have $i workers in the factory! Be careful not to break employee protection laws!")
    end

end