# Produce graph for the application section

using Statistics
using DelimitedFiles

include("generate_samples_Dirichlet.jl")
include("distance_Wasserstein.jl")
include("new_distance.jl")

# Which expriment to do
# 1. "N" means variations only in N, 
# 2. "alpha" means variations in alpha 
mode = "alpha"


# Law P_0: uniform over [0,1]
probability = rand 

# Upper bounds with the uniform over [0,1] as base measure 
# P1  
# 1/sqrt(N) * pi/8
# P_2
# (alphaV/(alphaV+1))^N * 1/3
# P_3
# 1/sqrt(N) * pi/8 * sqrt(alphaV/(alphaV+1))

# First we do only variations in N, and alpha is fixed 


if mode == "N"
    
    # alpha is fixed 
    global alphaV = 50. 

    # N is varying 
    NArray = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350, 400]
    nValues = length(NArray)


end 

if mode == "alpha"

    # alpha is varying  
    alphaArray = [10, 20, 30, 40 , 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
    nValues = length(alphaArray)

    # N is fixed 
    N = 50 

end



# Parameters for the approximation 
nTopNewDistance = 512 
nTopLowerBound = 5000
# Only used for the ground truth 
nBottomGT = 5000 
# Number of repetitions 
nRep = 10
# Number of point on the discrete grid 
nDiscrete = 200

# Build the output array 
outputArray = zeros(nValues,10)

for i = 1:nValues

    if mode == "N"
        global N = NArray[i]
        outputArray[i,1] = N 
        println("---------------")
        println("N = ", N)
        println("---------------")
    end 

    if mode == "alpha"
        global alphaV = alphaArray[i]
        println("---------------")
        println("alpha = ", alphaV)
        println("---------------")
        outputArray[i,1] = alphaV 
    end

    

    # Store the upper bounds 
    println("computing the upper bounds")
    outputArray[i,4] = 1/sqrt(N) * pi/8
    outputArray[i,7] = (alphaV/(alphaV+1))^N * 1/3
    outputArray[i,10] = 1/sqrt(N) * pi/8 * sqrt(alphaV/(alphaV+1))

    # Compute the lower bounds 
    println("computing the lower bounds")

    # nBottomGT = 5000
    
    gt = generateArrayDirichletProcess(nTopLowerBound,nBottomGT,alphaV,probability,"dependent")
    p1 = generateArrayDirichletProcess(nTopLowerBound,N,alphaV,probability,"finite")
    p2 = generateArrayDirichletProcess(nTopLowerBound,N,alphaV,probability,"stick")
    p3 = generateArrayDirichletProcess(nTopLowerBound,N,alphaV,probability,"dependent")

    outputArray[i,2] = lowerBound(gt, p1)
    outputArray[i,5] = lowerBound(gt, p2)
    outputArray[i,8] = lowerBound(gt, p3)

    # Approximate the new distance with nRep reptitions 
    println("computing the new distance")

    # nTopNewDistance=512
    auxArray = zeros(nRep,3)

    for j = 1:nRep 

        println("Repetition: ", j)

        gt = generateArrayDirichletProcess(nTopNewDistance,nBottomGT,alphaV,probability,"dependent")
        p1 = generateArrayDirichletProcess(nTopNewDistance,N,alphaV,probability,"finite")
        p2 = generateArrayDirichletProcess(nTopNewDistance,N,alphaV,probability,"stick")
        p3 = generateArrayDirichletProcess(nTopNewDistance,N,alphaV,probability,"dependent")

        auxArray[j,1] = dlip(gt, p1, 0. ,1. , nDiscrete)[1]
        auxArray[j,2] = dlip(gt, p2, 0. ,1. , nDiscrete)[1]
        auxArray[j,3] = dlip(gt, p3, 0. ,1. , nDiscrete)[1]

    end 

    outputArray[i,3] = Statistics.mean(auxArray[:,1])
    outputArray[i,6] = Statistics.mean(auxArray[:,2])
    outputArray[i,9] = Statistics.mean(auxArray[:,3])



end

# Store in the file ready to be read by pfgplot
# Only add manually 
# alpha mode 
# "alpha P1lb P1nd P1ub P2lb P2nd P2ub P3lb P3nd P3ub"
# N mode 
# "N P1lb P1nd P1ub P2lb P2nd P2ub P3lb P3nd P3ub"
# in the first line 

if mode == "N"
    writedlm("application_N.txt", outputArray)
end

if mode == "alpha"
    writedlm("application_alpha.txt", outputArray)
end