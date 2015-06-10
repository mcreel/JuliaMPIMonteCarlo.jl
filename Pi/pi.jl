reps = Int(1e6)  # desired number of MC reps
results = zeros(reps,1)
for sofar = 1:reps
    results[sofar,:] = 4.*(norm(rand(2,1)) .< 1.)
end    
println("pihat: ", mean(results))
