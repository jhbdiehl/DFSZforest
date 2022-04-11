


using CUDA
using BenchmarkTools
using SparseArrays
using CUDA.CUSPARSE


@btime cuA \ cux

e = 10

# Do it on CPU
A = [sparse(rand(Float64, 7, 7)) for i in 1:e]
x = rand(Float64, 7*e)

spA = blockdiag(A...)

@btime spA \ x

cuA = (CUDA.sparse(CUDA.rand(Float64, 7,7)), CUDA.sparse(CUDA.rand(Float64, 7,7)))


cuA = CUDA.rand(Float64, 7,7)
cux = CUDA.rand(Float64, 7)

@btime cuA \ cux

cu2A = [sparse(rand(Float64, 7, 7)) for i in 1:e]
cuspA = CUDA.blockdiag(A...)
# AHWEOAPRIHWOEIRHJNWOEFLWKJENF~~!~!!!!!!!!!!!!!!!!!!
blockdiag(CUDA.sparse(CUDA.rand(Float64, 7, 7)), CUDA.sparse(CUDA.rand(Float64, 7, 7)))

@btime cuspA \ cux


if Sys.iswindows()
    function hello_world()
        @cuprintf("Greetings from block %lld, thread %lld!\n", Int64(blockIdx().x), Int64(threadIdx().x))
        return
    end
else
    function hello_world()
       @cuprintf("Greetings from block %ld, thread %ld!\n", Int64(blockIdx().x), Int64(threadIdx().x))
       return
   end
end
@cuda blocks=2 threads=2 hello_world()
synchronize()



using CUDA: i32


function haversine_cpu(lat1::Float32, lon1::Float32, lat2::Float32, lon2::Float32, radius::Float32)
    c1 = cospi(lat1 / 180.0f0)
    c2 = cospi(lat2 / 180.0f0)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    d1 = sinpi(dlat / 360.0f0)
    d2 = sinpi(dlon / 360.0f0)
    t = d2 * d2 * c1 * c2
    a = d1 * d1 + t
    c = 2.0f0 * asin(min(1.0f0, sqrt(a)))
    return radius * c
end

function pairwise_dist_cpu(lat::Vector{Float32}, lon::Vector{Float32})
    # allocate
    n = length(lat)
    rowresult = Array{Float32}(undef, n, n)

    # brute force fill in each cell
    for i in 1:n, j in 1:n
        @inbounds rowresult[i, j] = haversine_cpu(lat[i], lon[i], lat[j], lon[j] , 6372.8f0)
    end

    return rowresult
end

# from https://devblogs.nvidia.com/parallelforall/fast-great-circle-distance-calculation-cuda-c/
function haversine_gpu(lat1::Float32, lon1::Float32, lat2::Float32, lon2::Float32, radius::Float32)
    # XXX: need to prefix math intrinsics with CUDA
    c1 = CUDA.cospi(lat1 / 180.0f0)
    c2 = CUDA.cospi(lat2 / 180.0f0)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    d1 = CUDA.sinpi(dlat / 360.0f0)
    d2 = CUDA.sinpi(dlon / 360.0f0)
    t = d2 * d2 * c1 * c2
    a = d1 * d1 + t
    c = 2.0f0 * CUDA.asin(CUDA.min(1.0f0, CUDA.sqrt(a)))
    return radius * c
end

# pairwise distance calculation kernel
function pairwise_dist_kernel(lat::CuDeviceVector{Float32}, lon::CuDeviceVector{Float32},
                              rowresult::CuDeviceMatrix{Float32}, n)
    i = (blockIdx().x-1i32) * blockDim().x + threadIdx().x
    j = (blockIdx().y-1i32) * blockDim().y + threadIdx().y

    if i <= n && j <= n
        # store to shared memory
        shmem = CuDynamicSharedArray(Float32, 2*blockDim().x + 2*blockDim().y)
        if threadIdx().y == 1
            shmem[threadIdx().x] = lat[i]
            shmem[blockDim().x + threadIdx().x] = lon[i]
        end
        if threadIdx().x == 1
            shmem[2*blockDim().x + threadIdx().y] = lat[j]
            shmem[2*blockDim().x + blockDim().y + threadIdx().y] = lon[j]
        end
        sync_threads()

        # load from shared memory
        lat_i = shmem[threadIdx().x]
        lon_i = shmem[blockDim().x + threadIdx().x]
        lat_j = shmem[2*blockDim().x + threadIdx().y]
        lon_j = shmem[2*blockDim().x + blockDim().y + threadIdx().y]

        @inbounds rowresult[i, j] = haversine_gpu(lat_i, lon_i, lat_j, lon_j, 6372.8f0)
    end

    return
end

function pairwise_dist_gpu(lat::Vector{Float32}, lon::Vector{Float32})
    # upload
    lat_gpu = CuArray(lat)
    lon_gpu = CuArray(lon)

    # allocate
    n = length(lat)
    rowresult_gpu = CuArray(zeros(Float32, n, n))

    # calculate a 2D block size from the suggested 1D configuration
    # NOTE: we want our launch configuration to be as square as possible,
    #       because that minimizes shared memory usage
    function get_threads(threads)
        threads_x = floor(Int, sqrt(threads))
        threads_y = threads ÷ threads_x
        return (threads_x, threads_y)
    end

    # calculate the amount of dynamic shared memory for a 2D block size
    get_shmem(threads) = 2 * sum(threads) * sizeof(Float32)

    kernel = @cuda launch=false pairwise_dist_kernel(lat_gpu, lon_gpu, rowresult_gpu, n)
    config = launch_configuration(kernel.fun, shmem=threads->get_shmem(get_threads(threads)))

    # convert to 2D block size and figure out appropriate grid size
    threads = get_threads(config.threads)
    blocks = ceil.(Int, n ./ threads)
    shmem = get_shmem(threads)

    kernel(lat_gpu, lon_gpu, rowresult_gpu, n; threads=threads, blocks=blocks, shmem=shmem)
    return Array(rowresult_gpu)
end

using Test

# generate reasonable data
function main(n = 1000)
    lat = rand(Float32, n) .* 45
    lon = rand(Float32, n) .* -120

    @test pairwise_dist_cpu(lat, lon) ≈ pairwise_dist_gpu(lat, lon) rtol=1e-2
end
main()