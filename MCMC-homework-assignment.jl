# Let's start by installing the Cosmology package!
using Pkg
Pkg.add("Cosmology")

# We'll also end up using all our old friends:
using WGLMakie
using CSV
using DataFrames
using Cosmology
using Statistics

# There is a data file in this directory, taken basically straight out of the Perlmutter+1999 paper.  We can read it with the CSV package.
data = CSV.read("p99-data.txt", DataFrame, delim=" ", ignorerepeated=true);

# Make a copy of the data columns that we want to treat as the "y" measurements.
# These are the measured brightnesses, and their Gaussian uncertainties (standard deviations).
data.mag = data.m_b_eff
data.sigma_mag = data.sigma_m_b_eff;

f = Figure()
Axis(f[1,1], title="Perlmutter+99 Supernovae", xlabel="Redshift z", ylabel="m_B")
errorbars!(data.z, data.mag, data.sigma_mag)
scatter!(data.z, data.mag, markersize=5, color=:maroon)
save("perlmutterPlot.png", f)
display(f)

# Here is how we will use the "cosmology" package.  This will create a cosmology "object" with the parameters we pass in.
# It does not take an Omega_Lambda parameter; instead, it takes Omega_Matter, and Omega_K (for "curvature"), where
# Omega_K = 1. - Omatter - Olambda.  We will also pass in "Tcmb=0", which tells it to ignore the effects of radiation.

universe = cosmology(OmegaK=0.1, OmegaM=0.4, Tcmb=0)
@show universe
@show universe.Ω_Λ;

# We can then pass that "universe" object to other functions to compute things about it.  Basically the only one you'll
# need is this `distance_modulus`, which tell you, in _magnitudes_, how much fainter an object is at the given redshift,
# versus how faint it would be if it were 10 parsecs away.

function distance_modulus(universe, z)
    DL = luminosity_dist(universe, z)
    # DL is in Megaparsecs; the distance for absolute to observed mag is 10 pc.
    5. * log10.(DL.val * 1e6 / 10.)
end;

# We'll cheat a bit and use a "nominal" cosmology with currently-accepted values of Omega_M = 0.29, Omega_DE = 0.71.
nominal = cosmology(Tcmb=0)

f = Figure()
ax = Axis(f[1,1], title="Perlmutter+99 Supernovae", xlabel="Redshift z", ylabel="Observed mag")
errorbars!(data.z, data.mag, data.sigma_mag)
scatter!(data.z, data.mag, markersize=5, color=:maroon)

# Compute the average absolute magnitude M given nominal cosmology -- ie, an estimate of the absolute mag of the supernovae
DLx = map(z->distance_modulus(nominal, z), data.z)
abs_mag = median(data.mag - DLx)

# Here's another way to plot a function evaluated on a grid of values.
zgrid = 0.01:0.01:1.
DL = map(z->distance_modulus(nominal, z), zgrid)
lines!(zgrid, DL .+ abs_mag, label="Nominal OmegaM = 0.29, OmegaDE = 0.71")

universe = cosmology(OmegaK=0.0, OmegaM=0.6, Tcmb=0)
DL = map(z->distance_modulus(universe, z), zgrid)
lines!(zgrid, DL .+ abs_mag, color=:red, label="OmegaM = 0.6, OmegaDE = 0.4")

universe = cosmology(OmegaK=0.0, OmegaM=0.1, Tcmb=0)
DL = map(z->distance_modulus(universe, z), zgrid)
lines!(zgrid, DL .+ abs_mag, color=:green, label="OmegaM = 0.1, OmegaDE = 0.9")

#f[2,1] = Legend(f, ax, "Cosmologies", framevisible = false)
# Create a legend for our plot
axislegend(ax, position = :rb)
save("PerlmutterLines.png", f)
f

# Here's our scalar estimate of the absolute mag.
abs_mag

function supernova_log_likelihood(zs, mag, mag_error, M, Omatter, Ode)
    # z: vector of redshifts
    # mag: vector of measured magnitudes
    # mag_error: vector of uncertainties on the measured magnitudes (sigmas).
    # M: scalar, absolute magnitude of a Type-1a supernova
    # Omatter: scalar Omega_M, amount of matter in the universe
    # Ode: scalar Omega_DE, amount of dark energy in the universe

    ###   YOUR CODE HERE!!
    cosModel = cosmology(OmegaK=1. - Omatter - Ode, OmegaM=Omatter, Tcmb=0) #make the cosmology "object"

    mag_pred = (z -> distance_modulus(cosModel, z) + M).(zs) #vectorize the operations on zs

    chi = (mag_pred - mag) ./ mag_error
    sum(-0.5 .* chi.^2)

    # You must return a scalar value
end;

n_om, n_ode = 50,50
om_grid = LinRange(0, 1, n_om)
ode_grid = LinRange(0, 1, n_ode)
sn_ll = zeros(n_om, n_ode)
for i in 1:n_om
    for j in 1:n_ode
        sn_ll[i, j] = supernova_log_likelihood(data.z, data.mag, data.sigma_mag, abs_mag, om_grid[i], ode_grid[j])
    end
end
f = Figure()
ax = Axis(f[1, 1])

heatmap!(om_grid, ode_grid, sn_ll, colorrange=[maximum(sn_ll)-20, maximum(sn_ll)])

save("FirstHeatmapLimitedRange.png", f)
f

n_om2, n_ode2 = 100,100
om_grid2 = LinRange(0, 3, n_om2)
ode_grid2 = LinRange(-1.5, 3, n_ode2)
sn_ll2 = zeros(n_om2, n_ode2)
for i in 1:n_om2
    for j in 1:n_ode2
        try #very bad way of doing this, should rather make supernova_log_likelihood to check for parameters.
            sn_ll2[i, j] = supernova_log_likelihood(data.z, data.mag, data.sigma_mag, abs_mag, om_grid2[i], ode_grid2[j])
        catch err
            sn_ll2[i, j] = -Inf
        end
    end
end

f = Figure()
ax = Axis(f[1, 1])

heatmap!(om_grid2, ode_grid2, sn_ll2, colorrange=[maximum(sn_ll2)-20, maximum(sn_ll2)])

save("HeatmapExtendedRange.png", f)
f

function cornerplot(x, names; figsize=(600,600))
    # how many columns of data
    dim = size(x, 2)
    # rows to plot
    idxs = 1:size(x,1)
    f = Figure(size=figsize)
    for i in 1:dim, j in 1:dim
        if i < j
            continue
        end
        ax = Axis(f[i, j], aspect = 1,
                  topspinevisible = false,
                  rightspinevisible = false,)
        if i == j
            hist!(x[idxs,i], direction=:y)
            ax.xlabel = names[i]
        else
            #scatter!(x[idxs,j], x[idxs,i], markersize=4)
            hexbin!(x[idxs,j], x[idxs,i])
            ax.xlabel = names[j]
            ax.ylabel = names[i]
        end
    end
    f
end;

function mcmc_cyclic(logprob_func, propose_func, initial_p, n_steps)
    p = initial_p
    logprob = logprob_func(p)
    chain = zeros(n_steps, length(p))
    n_accept = zeros(length(p))

    for i in 1:n_steps

        # We're going to update one index at a time... 1, 2, 1, 2, ....
        update_index = 1 + ((i-1) % length(p))

        # Call the proposal function to generate new values for all parameters...
        p_prop = propose_func(p)
        # ... but then only keep one of the new parameter values!
        p_new = copy(p)
        p_new[update_index] = p_prop[update_index]
        
        logprob_new = logprob_func(p_new)

        ratio = exp(logprob_new - logprob)
        if ratio > 1
            # Jump to the new place
            p = p_new
            logprob = logprob_new
            n_accept[update_index] += 1
        else
            # Jump to the new place with probability "ratio"
            u = rand()
            if u < ratio
                # Jump to the new place
                p = p_new
                logprob = logprob_new
                n_accept[update_index] += 1
            else
                # Stay where we are
            end
        end
        chain[i, 1:end] = p
    end
    # The number of times we step each parameter is roughly n_steps / n_parameters
    chain, n_accept ./ (n_steps ./ length(p))
end;

function propose(p, jump_sizes)
    p .+ randn(length(p)) .* jump_sizes
end;

initial_guess = [abs_mag, .5, .5]
jump_sizes = [0.02, 0.09, 0.1]

chain, accept_rate = mcmc_cyclic(p -> 
        supernova_log_likelihood(data.z, data.mag, data.sigma_mag, p[1], p[2], p[3]),
p -> propose(p, jump_sizes),
initial_guess, 150000)
println("Acceptance Rates")
accept_rate

cplot1 = cornerplot(chain, ["M", "Omega_M", "Omega_DE"])
save("cornerplotFirstMCMC.png", cplot1)

f = Figure(size = (1500, 400))
ax = Axis(f[1,1], title="MCMC chain for M", xlabel="step", ylabel="M")

scatter!(chain[1:10:end, 1])
ax = Axis(f[1,2], title="MCMC chain for Omega_M", xlabel="step", ylabel="Omega_M")
scatter!(chain[1:10:end, 2])
ax = Axis(f[1,3], title="MCMC chain for Omega_DE", xlabel="step", ylabel="Omega_DE")
scatter!(chain[1:10:end, 3])

save("MCMCChain.png", f)
display(f)

OmegaMMCMC = chain[begin:end, 2]
OmegaDEMCMC = chain[begin:end, 3];

n_ombins, n_odebins = 20,20

om_gridMCMC = LinRange(minimum(OmegaMMCMC), maximum(OmegaMMCMC), n_ombins + 1)
ode_gridMCMC = LinRange(minimum(OmegaDEMCMC), maximum(OmegaDEMCMC), n_odebins + 1)

minOmegaMMCMC = minimum(OmegaMMCMC)
minOmegaDEMCMC = minimum(OmegaDEMCMC)
maxOmegaMMCMC = maximum(OmegaMMCMC)
maxOmegaDEMCMC = maximum(OmegaDEMCMC)
rangeOmegaMMCMC = maxOmegaMMCMC - minOmegaMMCMC
rangeOmegaDEMCMC = maxOmegaDEMCMC - minOmegaDEMCMC

histogramValues = zeros(n_ombins + 1, n_odebins + 1)
for i in 1:size(OmegaDEMCMC)[begin]
    histogramValues[trunc(Int, (OmegaMMCMC[i] - minOmegaMMCMC) * n_ombins / (rangeOmegaMMCMC)) + 1, trunc(Int, (OmegaDEMCMC[i] - minOmegaDEMCMC) * n_odebins / rangeOmegaDEMCMC) + 1] += 1
## I was doing to more manually, which is probably correct, but very slow.
#    for j in 1:(n_ombins - 1)
#        for k in 1:(n_odebins - 1)
#            histogramValues[j, k] += (om_gridMCMC[j] <= OmegaMMCMC[i] && OmegaMMCMC[i] <= om_gridMCMC[j + 1] && ode_gridMCMC[k] <= OmegaDEMCMC[i] && OmegaDEMCMC[i] <= ode_gridMCMC[k + 1])
#        end
#    end
end
histogramValues;

f = Figure()
Axis(f[1, 1])

contour!(om_gridMCMC, ode_gridMCMC, histogramValues)

save("MyHistogram.png", f)

display(f)

Pkg.add(["FHist", "AffineInvariantMCMC"])


using FHist

h = FHist.Hist2D((chain[:,2], chain[:,3]); nbins=(100,100))

counts = bincounts(h);
xc,yc = bincenters(h);

f = Figure()
Axis(f[1, 1])

contour!(xc, yc, counts) #, levels=[10,50,100])
#really small/many bins make this contour plot really jaggered

save("MyFirstJaggeredContour.png", f)
f

h = FHist.Hist2D((chain[:,2], chain[:,3]); nbins=(20,20))

counts = bincounts(h);
xc,yc = bincenters(h);

f = Figure()
Axis(f[1, 1])

contour!(xc, yc, counts)#, levels=[10,50,100])

save("SmallerBinsCountour.png", f)
f

using AffineInvariantMCMC

numdims = 3
numwalkers = 50
thinning = 10
numsamples_perwalker = 15000
burnin = 5000

# Start out by doing a "burn-in" short run...
initial = [abs_mag, 0.5, 0.5] .+ randn(numdims, numwalkers)*0.01
chain, ll = AffineInvariantMCMC.sample( p -> 
        try
            supernova_log_likelihood(data.z, data.mag, data.sigma_mag, p[1], p[2], p[3])
        catch err
            -Inf
        end, numwalkers, initial, burnin, 1)

# Then start the "main" run where the burn-in finished:
initial = chain[:,:,end]
chain, ll = AffineInvariantMCMC.sample( 
   p -> 
        try
            supernova_log_likelihood(data.z, data.mag, data.sigma_mag, p[1], p[2], p[3])
        catch err
            -Inf
        end
    , numwalkers, initial, numsamples_perwalker, thinning)
flatchain, flatllhoodvals = AffineInvariantMCMC. flattenmcmcarray(chain, ll);

# This "flatchain" is transposed from the way we produced it in our code, so you can do

chain = flatchain';

cplot2 = cornerplot(chain, ["M", "Omega_M", "Omega_DE"])

save("AffineMCMCCPlot.png", cplot2)
cplot2

h = FHist.Hist2D((chain[:,2], chain[:,3]); nbins=(20,20))

counts = bincounts(h);
xc,yc = bincenters(h);

f = Figure()
Axis(f[1, 1])

contour!(xc, yc, counts)#, levels=[10,50,100])

save("AffineMCMCContour.png", f)
f

h = FHist.Hist2D((chain[:,2], chain[:,3]); nbins=(100,100))

counts = bincounts(h);
xc,yc = bincenters(h);

f = Figure()
Axis(f[1, 1])

contour!(xc, yc, counts) #, levels=[10,50,100])

save("finalsmoothcountour.png", f)
f
