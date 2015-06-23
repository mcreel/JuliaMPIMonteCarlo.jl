include("Auction.jl") # load auction model code
include("AIS.jl") # the adaptive importance sampling algorithm

function AuctionWrapper()
    # true theta
    theta = [0.5 0.5]
    # generate 'true' aux. stat.
    Zn = aux_stat(theta)
    nParticles = 500 # particles to keep per iter
    multiples = 5  # particles tried is this multiple of particle kept
    StopCriterion = 0.1 # stop when proportion of new particles accepted is below this
    AISdraws = 5000
    neighbors = 25
    contrib = AIS_fit(Zn, nParticles, multiples, StopCriterion, AISdraws, neighbors)
    println("ABC estimator:", contrib)
end

for i = 1:1
    AuctionWrapper()
end    
