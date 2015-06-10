# The auction model, inputs are parameters and sample size
function AuctionModel(theta::Array{Float64,2}, n::Int64, setseed=false)
	
    # SMIL requires same random draws
    if setseed
        srand(1234)
    end    
    
    # the model
	theta1 = theta[:,1]
	theta2 = theta[:,2]
	N = 6
	# quality of good
	x = rand(n,1)
	# valuations drawn from exponetial mean phi
	phi = exp(theta1 .+ theta2.*x)
	# highest valuation
	v = -log(minimum(rand(n,6),2)).*phi
	# get winning bid
	z = v./phi
	D = exp(-5.*z).*(60.*exp(5.*z) + 300.*phi .* exp(4.*z) - 300.*phi .* exp(3.*z)
	    + 200.*phi .* exp(2.*z) - 75.*phi .* exp(z) + 12.*phi)/60. - 137.*phi/60.
	b = v - D ./ ((1. - exp(-v./phi)).^(N-1))
	b = b.*(b.>0)
	data = [b x]
end

# Auxiliary statistic.
function aux_stat(theta::Array{Float64,2})
    n = 80
	data = AuctionModel(theta, n)
	b = data[:,1]
	# bound bid for numeric stability
	b = (b.>0.01).*b + (b.<0.01)*0.01
	y = log(b) # use log of bid
	n = size(y,1)
    x = [ones(n,1) data[:,2]]
	bhat = x\y
	e = y - x*bhat
	sig = log(e'*e/(n-2.))
	m1 = mean(log(b))
	m2 = std(log(b))
	m3 = mean((log(b)-m1).^3)
	Z = [bhat' sig m1 m2 m3]
    return Z
end

# this function generates a draw from the prior
function sample_from_prior()
	theta = rand(1,2)
    lb = [-1 0]
    ub = [3 2]
    theta = (ub-lb).*theta + lb
end

# the prior: needed to compute AIS density, which uses
# the prior as a mixture component, to maintain the support
function prior(theta::Array{Float64,2})
    lb = [-1 0]
    ub = [3 2]
    c = 1./prod(ub - lb)
    p = ones(size(theta,1),1)*c
    ok = all((theta.>=lb) & (theta .<=ub),2)
    p = p.*ok
end    

function check_in_support(theta::Array{Float64,2})
    lb = [-1 0]
    ub = [3 2]
    ok = (theta .>= lb) & (theta .<= ub)
    ok = mean(ok)==1.
    return ok, lb, ub
end

