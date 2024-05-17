# Produce the convergence graph for Wasserstein over Wasserstein  
# or for our new distance 

using Statistics
using DelimitedFiles

include("generate_samples_Dirichlet.jl")
include("distance_Wasserstein.jl")
include("new_distance.jl")

# -------------------------------------------------------------------
# Parameters 
# -------------------------------------------------------------------

# Decide which distance to compute among  
# 1. "WoW" Wasserstein over Wasserstein [To use only with the dependent Dirichlet process] 
# 2. "new_distance" the HIPM distance dlip 
# 3. "lower_bound" the lower bound given by (1)
toCompute = "new_distance"

# Decide which input measures to use 
# 1. "same" (bottom of the figure)
# 2. "splitting" (top of the figure)
baseMeasure = "same"

# Exponent for the Wasserstein distance 
p = 1.

# Concentration parameter for the Dirichlet process 
alpha = 1. 

# Parameter for the new HIPM distance 
if baseMeasure == "same"
    a,b = 0. , 1.
end 
if baseMeasure == "splitting"
    a,b = -1., 1.
end

nDiscrete = 150


# --------------------------------------------------------------------
# Define base measures 
# --------------------------------------------------------------------


function probability1()

    if baseMeasure == "same"
        return rand()
    elseif baseMeasure == "splitting"
        return rand() .- 0.5
    end

end

function probability2()

    if baseMeasure == "same"
        return rand()
    elseif baseMeasure == "splitting"
        atom = rand()
        mixture = rand((0,1))
        return mixture * ( -1. + 0.25 * atom ) + (1 - mixture) * (0.75 + 0.25 * atom)
    end

end

# -------------------------------------------------------------------
# Function to do the stats  
# -------------------------------------------------------------------

function empiricalDistAverage(alpha, nTop, nBottom,nRep,toCompute)
    # nTop is the outer number 
    # nBottom is the inner number 
    # nRep is how many time we repeat to do the stats 

    arrayAnswer = zeros(nRep)

    # Threads.@threads for i in 1:nRep
    for i in 1:nRep 

        println("Computation of one distance")

        # Generate the two arrays of DP
        measure1 = generateArrayDirichletProcess(nTop,nBottom,alpha,probability1,"dependent")
        measure2 = generateArrayDirichletProcess(nTop,nBottom,alpha,probability2,"dependent")


        # Compute the answer
        if toCompute == "WoW"
            # To use only for the dependent Dirichlet process
            arrayAnswer[i] = wassersteinOverWasserstein(measure1[:,:,1], measure2[:,:,1], p)
        end 

        if toCompute == "new_distance"
            arrayAnswer[i] = dlip(measure1, measure2, a,b, nDiscrete,1000,5,1e-4)[1]
        end

        if toCompute == "lower_bound"
            arrayAnswer[i] = lowerBound(measure1,measure2)
        end 

    end

    return Statistics.mean(arrayAnswer), Statistics.std(arrayAnswer)

end

# -------------------------------------------------------------------
# Parameters to do the stats 
# -------------------------------------------------------------------

nTopArray = [16,32,64,128,256]
nBottom = 5000
nRep = 24

nArray = length(nTopArray)
outputArray = zeros(nArray,3)

for i = 1:nArray

    nTop = nTopArray[i]

    println("-"^40)
    println("nTop ", nTop)
    println("-"^40)

    outputArray[i,1] = nTop
    outputArray[i,2],outputArray[i,3] = empiricalDistAverage(alpha, nTop, nBottom,nRep,toCompute)
    
end 

println(outputArray)

# Store in the file ready to be read by pfgplot
# Only add manually "n mean error" in the first line
writedlm("numerics_"*toCompute*"_"*baseMeasure*".txt", outputArray)