# WHEN DOING AIS fit, there is no check that the theta sampled from AIS density
# is in the support of the prior. This is not a problem, unless the support
# restriction is needed to avoid crashes. Also, it's inefficient to
# compute the statistic if the parameter draw is out of the support
# FIX THIS, eventually

using Distributions
using Distances
#=
AIS. The adaptive method is used to find the particles that
will be used to construct the importance sampling density.
The particles are found using perturbations of single 
parameters, and the particles are kept withing user provided
bounds. However, the importance sampling density that uses
the particles is a mixture of normals, each centered on a
particle, and the support is R^k. Thus, no truncation 
occurs in the importance sampling density.
=# 

# sample from particles: normal perturbation of a 
# single parameter, with bounds enforcement by rejection
# sampling
function sample_from_particles(particles::Array{Float64,2}, delta::Array{Float64,2} )
	n, k = size(particles)
    i = rand(1:n)
	j = rand(1:k)
    ok = 0.
    theta_s = similar(particles[i,:])
    @inbounds while ok != 1.
        theta_s = particles[i,:]
        theta_s[:,j] = theta_s[:,j] + delta[:,j]*randn()
        ok, junk, junk = check_in_support(theta_s)
    end    
    return theta_s
end

# sample from AIS density: choose random particle,
# add a MVN perturbation
function sample_from_AIS(particles::Array{Float64,2})
    delta = 0.1*std(particles,1)
	i = rand(1:size(particles,1))
    theta_s = particles[i,:]
    theta_s = theta_s + delta.*randn(size(delta))
    return theta_s
end

# the importance sampling density: mixture of normals
function AIS_density(theta::Array{Float64,2}, particles::Array{Float64,2})
    # what scaling to use here?
    delta  = 0.1*std(particles,1)
    sig = diagm(vec(delta.^2))
    nthetas = size(theta,1)
    nparticles = size(particles,1)
    dens = zeros(nthetas,1)
    @inbounds for i = 1:nparticles
        @inbounds for j = 1:nthetas
            thetaj = theta[j,:]
            thetaj2 = vec(thetaj)
            mu = vec(particles[i,:])
            d = MvNormal(mu, sig)
            dens[j,1] += pdf(d, thetaj2)
        end
    end
    dens = dens/nparticles
    return dens
end    

function AIS_fit(Zn, nParticles, multiples, StopCriterion, AISdraws, nn="default", whichfit="ll")
    # neighbors
    if nn=="default"
        neighbors = round(Int, 5.0*AISdraws^0.25)
    else
        neighbors = nn
    end  
    # do AIS to get particles
    particles, Zs = AIS_algorithm(nParticles, multiples, StopCriterion, Zn)
    # sample from AIS particles
    thetas = zeros(AISdraws, size(particles,2))
    Zs = zeros(AISdraws, size(Zn,2))
    @inbounds for i = 1:AISdraws
        thetas[i,:] = sample_from_AIS(particles)
        Zs[i,:] = aux_stat(thetas[i,:])
    end    
    # compute scaling limiting outliers
    dimZ = size(Zs,2)
    Z2 = copy(Zs)
    @inbounds for i = 1:dimZ
        q = quantile(Z2[:,i],0.99)
        # top bound
        test =  Z2[:,i] .< q
        Z2[:,i] = Z2[:,i].*test  + q.*(1. - test)
        q = -quantile(-Z2[:,i],0.99)
        # bottom bound
        test =  Z2[:,i] .> q
        Z2[:,i] = Z2[:,i] .* test + q.*(1. - test)
    end
    stdZ = std(Z2,1)
    distances = pairwise(Euclidean(),(Zs./stdZ)', (Zn./stdZ)') # get all distances
    ind = sortperm(vec(distances)) # indices of k nearest neighbors
    selected = ind[1:neighbors]
    thetas = thetas[selected,:] # the nearest neighbors
    distances = distances[selected,:]
    Zs = Zs[selected,:]
    # do the AIS weighting, and drop out-of-support draws
    AISweights = prior(thetas) ./ AIS_density(thetas,particles)
    test = (AISweights .> 0.)
    thetas = thetas[test[:,1],:] # the nearest neighbors
    distances = distances[test[:,1],:]
    Zs = Zs[test[:,1],:]
    AISweights = AISweights[test[:,1],:] 
    # overall weights are AIS times kernel weight
    m = maximum(distances)
    if m > 0
        weight = 2.*distances/m
    else
        weight = 0.
    end    
    weight = AISweights.*pdf(Normal(),weight)
    weight = weight/sum(weight)
    # compute estimator: default is local linear, but local constant or both are available
    if whichfit=="lc"
        fit = sum(thetas.*weight,1)
        return fit
    elseif whichfit=="ll"
        X = [ones(size(Zs,1),1) Zs]
        XX = weight .* X;
        b = inv(X'*XX)*XX'*thetas
        fit = [1. Zn]*b
        return fit
    else    
        lc_fit = sum(thetas.*weight,1)
        X = [ones(size(Zs,1),1) Zs]
        XX = weight .* X;
        b = inv(X'*XX)*XX'*thetas
        ll_fit = [1. Zn]*b
        return [lc_fit ll_fit]
    end   
end    


# Draw particles from prior
function GetInitialParticles(nParticles::Int64)
    particle = sample_from_prior()
    Z = aux_stat(particle)
    Zs = zeros(nParticles, size(Z,2))
    particles = zeros(nParticles, size(particle,2))
    particles[1,:] = particle
    Zs[1,:] = Z
    @inbounds for i = 2:nParticles
        particles[i,:] = sample_from_prior()
        Z = aux_stat(particles[i,:])
        Zs[i,:]	= Z
    end
    return particles, Zs
end

# Sample from current particles
function GetNewParticles(particles::Array{Float64,2}, multiple::Int64)
    nParticles, dimTheta = size(particles)
    newparticles = zeros(multiple*nParticles,dimTheta)
    # get size of Z
    particle = sample_from_prior()
    Z = aux_stat(particle)
    dimZ = size(Z,2)
    newZs = zeros(multiple*nParticles, dimZ)
    delta = 0.5*std(particles,1)
    @inbounds for i = 1:multiple*nParticles
        newparticles[i,:] = sample_from_particles(particles, delta)
        newZs[i,:] = aux_stat(newparticles[i,:])
    end
    return newparticles, newZs
end

# Select the best particles from current and new
function Select(nParticles, distances::Array{Float64,1}, particles::Array{Float64,2}, Zs::Array{Float64,2})
    ind = sortperm(distances) # indices of distances
    ind = ind[1:nParticles] # indices of best
    new = sum(ind .> nParticles)/nParticles # percentage of new particles kept
    # keep the best distances, particles, Zs
    distances = distances[ind]
    particles = particles[ind,:]
    Zs = Zs[ind,:]
    return new, particles, Zs, distances
end    

#=
function GetPlot(iter::Int64, particles::Array{Float64,2}, axes, ax)        
    if iter == 1
        fig, axes = subplots(5, 3,sharex="all", sharey="all")
    end
    ax = axes[iter, 1]
    ax[:scatter](particles[:,1],particles[:,2])
    ax = axes[iter, 2]
    ax[:scatter](particles[:,1],particles[:,3])
    ax = axes[iter, 3]
    ax[:scatter](particles[:,2],particles[:,3])
    return axes, ax
end        
=#

function AIS_algorithm(nParticles::Int64, multiple::Int64, StopCriterion::Float64, Zn::Array{Float64,2}, DoPlot::Bool=false)
    # the initial particles
    particles, Zs = GetInitialParticles(nParticles)
    # do bounding to compute scale the first time
    # in the loop, selection will remove outliers
    dimZ = size(Zs,2)
    @inbounds for i = 1:dimZ
        q = quantile(Zs[:,i],0.99)
        # top bound
        test =  Zs[:,i] .< q
        Zs[:,i] = Zs[:,i].*test  + q.*(1. - test)
        q = -quantile(-Zs[:,i],0.99)
        # bottom bound
        test =  Zs[:,i] .> q
        Zs[:,i] = Zs[:,i] .* test + q.*(1. - test)
    end
    iter = 0
    # the main loop
    new = 1.
    @inbounds while mean(new) > StopCriterion # currently, stops when number of accepted new particles is low enough 
        iter +=1
        stdZ = std(Zs,1)
        # generate new particles
        newparticles, newZs = GetNewParticles(particles, multiple)
        # stack them
        Zs = [Zs; newZs]
        particles = [particles; newparticles]
        # find distances
        distances = pairwise(Euclidean(),(Zs./stdZ)', (Zn./stdZ)') # get all distances
        distances = vec(distances)
        # select from new and old
        new, particles, Zs, distances =  Select(nParticles, distances, particles, Zs)
        #=if DoPlot & (iter <=5)
            GetPlot(iter, particles, axes, ax)
        end
        =#
    end
    return particles, Zs
end
