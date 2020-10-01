# JuliaMPIMonteCarlo
Illustrates using Julia and MPI to do Monte Carlo. Release v1.0 contains the versions used in the paper "A Note on Julia and MPI, with Code Examples" (see http://link.springer.com/article/10.1007/s10614-015-9516-5). That code will no longer run.

# current info:
The code is currently being adapted to run on Julia 1.x. montecarlo.jl will be added back in, and the pi example resurrected. The auction code will be adapted, too.

# out of date info:
(this is no longer true) Since it was written, montecarlo.jl and the pi example have been incorporated in the MPI.jl package (https://github.com/JuliaParallel/MPI.jl), and the interface has changed slightly. The pi example has been deleted from this repo. The code for the Auction model adjusted to use the new montecarlo.jl interface is in the repo https://github.com/mcreel/ABCAuction.jl

