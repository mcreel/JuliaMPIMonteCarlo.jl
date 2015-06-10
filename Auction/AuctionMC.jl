include("Auction.jl") # load auction model code
include("AIS.jl") # the adaptive importance sampling algorithm
include("montecarlo.jl")

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
end

# the monitoring function
function AuctionMonitor(sofar, results)
    if mod(sofar,100) == 0
        theta = [0.5 0.5]
        m = mean(results[1:sofar,1:end],1)
        er = theta - m; # theta defined at top level, so ok to use
        b = mean(er,1)
        s = std(results[1:sofar,:],1) 
        mse = s.^2 + b.^2
        rmse = sqrt(mse)
        println()
        println("reps so far: ", sofar)
        println("mean: ", m)
        println("bias: ", b)
        println("st. dev.: ", s)
        println("rmse.: ",rmse)
    end
end

function main()
    reps = 100   # desired number of MC reps
    n_returns = 2
    pooled = 1  # do this many reps b
    montecarlo(AuctionWrapper, AuctionMonitor, reps, n_returns, pooled)
end

main()
