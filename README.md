# JuliaMPIMonteCarlo
Illustrates using Julia and MPI to do Monte Carlo. Release v1.0 contains the versions used in the paper "A Note on Julia and MPI, with Code Examples" (see http://link.springer.com/article/10.1007/s10614-015-9516-5). That code will no longer run. Since it was written, montecarlo.jl and the pi example have been incorporated in the MPI.jl package (https://github.com/JuliaParallel/MPI.jl), and the interface has changed slightly. The pi example has been deleted from this repo. The code for the Auction model adjusted to use the new montecarlo.jl interface is in the repo https://github.com/mcreel/ABCAuction.jl

Thus, this repo now just archives the code used in the paper, but for working code to do the same things, you should go to the sources mentioned above.
