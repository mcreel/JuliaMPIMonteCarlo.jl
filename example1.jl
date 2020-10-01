using Pkg
Pkg.activate("./")
using MPI, Statistics
function main()
    if !MPI.Initialized()
        MPI.Init()
    end    
    comm = MPI.COMM_WORLD
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    recv_mesg = zeros(100,1)
    if rank == 0
        for i = 1:size-1
            MPI.Send(rand(100), i, 0, comm)
        end
    else
        MPI.Recv!(recv_mesg, 0, 0, comm)
    end    
    MPI.Barrier(comm)
    m = mean(recv_mesg,dims=1)
    for i=1:size-1
        if rank==i
            println("Results on rank ",i,": ",m)
        end
        sleep(rank) # print out results in rank order
    end    
    MPI.Finalize()
end
main()
