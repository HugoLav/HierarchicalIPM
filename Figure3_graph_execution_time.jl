# Computing the execution times for different siutation


using Statistics
using DelimitedFiles

include("generate_samples_Dirichlet.jl")
include("distance_Wasserstein.jl")
include("new_distance.jl")

# Define the parameters 
# Three types of variations: in nTop, nBottom and nGrid 
# 1. "nTop", when nTop (n in the article) varies 
# 2. "nBottom", when nBottom (m in the article) varies
# 3. "nGrid", when nGrid (M in the article) varies
mode = "nTop"

alpha = 1
nRep = 12



# -------------------------------------------
# Define the values of parameters tested 
# -------------------------------------------

# Initialize the variables to avoid having local/global problem
nTop = 10
nBottom = 10
nGrid = 10

if mode == "nTop"

    nTopArray = [16,32,64,128,256]
    global nBottom = 5000 
    global nGrid = 250

    lengthArray = length(nTopArray)

end

if mode == "nBottom"

    global nTop = 128
    nBottomArray = [50,100,500,1000,2500,5000]
    global nGrid = 250

    lengthArray = length(nBottomArray)

end

if mode == "nGrid"

    global nTop = 128
    global nBottom = 1000
    nGridArray = [32,64,128,256,512]

    lengthArray = length(nGridArray)

end


# -------------------------------------------
# Create the array to store the output  
# -------------------------------------------

if mode == "nBottom" || mode == "nTop"

    outputArray = zeros(lengthArray,5)

end

if mode == "nGrid"

    outputArray = zeros(lengthArray,3)

end

# -------------------------------------------
# Function for the initial conditions   
# -------------------------------------------

# Same measures  

function probability1Same()
    return rand()
end

function probability2Same()
    return rand()
end

# Splitting example 

function probability1Split()
    return rand() .- 0.5
end

function probability2Split()
    atom = rand()
    mixture = rand((0,1))
    return mixture * ( -1. + 0.25 * atom ) + (1 - mixture) * (0.75 + 0.25 * atom)
end


# -------------------------------------------
# Run the loop    
# -------------------------------------------

# Do first some random computations to get the system to compile 

for k = 1:nRep 

    measure1 = generateArrayDirichletProcess(10,10,alpha,probability1Same,"dependent")
    measure2 = generateArrayDirichletProcess(10,10,alpha,probability2Same,"dependent")
        
    dlip(measure1, measure2, 0.,1., nGrid,1000,5,1e-4)
    wassersteinOverWasserstein(measure1[:,:,1], measure2[:,:,1], 1.)

end 

# Actual timing of the code 

for i = 1:lengthArray

    # Collect the parameter which varies 
    if mode == "nTop"
        global nTop = nTopArray[i]
    end 

    if mode == "nBottom"
        global nBottom = nBottomArray[i]
    end 

    if mode == "nGrid"
        global nGrid = nGridArray[i]
    end

    println("-----------------")
    println("nTop ", nTop)
    println("nBottom ", nBottom)
    println("nGrid ", nGrid)

    # Initialize the total time spent  
    sSameND = 0.
    sSplitND = 0. 
    sSameWoW = 0.
    sSplitWoW = 0. 

    for j = 1:nRep 

        # Initialize the measures: same 
        measure1 = generateArrayDirichletProcess(nTop,nBottom,alpha,probability1Same,"dependent")
        measure2 = generateArrayDirichletProcess(nTop,nBottom,alpha,probability2Same,"dependent")
        
        # Compute the distance 
        sSameND += @elapsed dlip(measure1, measure2, 0.,1., nGrid,1000,5,1e-4)
        if mode != "nGrid"
            sSameWoW += @elapsed wassersteinOverWasserstein(measure1[:,:,1], measure2[:,:,1], 1.)
        end 

        # Initialize the measures: split situation 
        measure1 = generateArrayDirichletProcess(nTop,nBottom,alpha,probability1Split,"dependent")
        measure2 = generateArrayDirichletProcess(nTop,nBottom,alpha,probability2Split,"dependent")
        
        # Compute the distance 
        sSplitND += @elapsed dlip(measure1, measure2, -1.,1., nGrid,1000,5,1e-4)
        if mode != "nGrid"
            sSplitWoW += @elapsed wassersteinOverWasserstein(measure1[:,:,1], measure2[:,:,1], 1.)
        end 

    end

    # Store the results 

    # Parameters which varies 
    if mode == "nTop"
        outputArray[i,1] = nTop
    end

    if mode == "nBottom"
        outputArray[i,1] = nBottom
    end 

    if mode == "nGrid"
        outputArray[i,1] = nGrid 
    end 

    # Then store the results for the new distance 
    outputArray[i,2] = sSameND / nRep
    outputArray[i,3] = sSplitND / nRep 
    
    # If we are not on nGrid store for WoW 
    if mode != "nGrid"
        outputArray[i,4] = sSameWoW / nRep
        outputArray[i,5] = sSplitWoW / nRep 
    end


end

println(outputArray)


# Save the output in a txt file ready to be read by pfgplot
# What to copy paste on top of the txt file 
# 1. nTop 
# nTop sameND splitND sameWoW splitWoW
# 2. nBottom 
# nBottom sameND splitND sameWoW splitWoW
# 3. nGrid 
# nGrid sameND splitND 

writedlm("timing_"*mode*".txt", outputArray)