# Generating samples from a Dirichlet process
# Include all three approximation of Section 6. 
# The main function is 
# generateArrayDirichletProcess
# which generates realization of \tilde{\mathhb{Q}}_{(n,m)} for the tree approximations.  


using Distributions


# Measures are given as a nTop x nBottom x 2 array
# nTop = n in the article 
# nBottom = m in the article 
# nGrid = M in the article 
# a[:,:,1] -> location of the atom 
# a[:,:,2] -> mass of the atom 




function generateDirichletProcessWithoutWeight(N,alpha,probability)
    # N number of samples 
    # alpha concentration parameter 
    # P is the function to generate the samples 
    # return list of atoms 

    if N == 0
        return []
    else
        previous = generateDirichletProcessWithoutWeight(N-1,alpha,probability)
        if rand() <= alpha/(alpha+N-1)
            return push!(previous, probability())
        else
            index = rand(1:(N-1))
            return push!(previous, previous[index])
        end
    end  

end

function generateDependentDirichletProcess(N,alpha,probability)

    output = zeros(N,2)
    output[:,1] = generateDirichletProcessWithoutWeight(N,alpha,probability)
    output[:,2] = fill(1/N,N)

    return output

end 

function generateStickBreaking(N,alpha,probability)

    # First stick breaking of the Dirichlet process 
    output = zeros(N,2)
    stick = 1. 

    # Do the loop 
    for i = 1:(N-1) 
        output[i,1] = probability()
        weight = rand(Beta(1., alpha))
        output[i,2] = weight * stick 
        stick -= weight * stick 
    end

    output[end,1] = probability()
    output[end,2] = stick 

    return output


end

function generateFiniteDirichletProcess(N,alpha,probability)

    output = zeros(N,2)

    for i=1:N
        output[i,1] = probability()
    end

    output[:,2] = rand( Dirichlet(N, alpha/N) )

    return output 

end

function generateArrayDirichletProcess(nTop,nBottom,alpha,probability,mode)

    output = zeros(nTop,nBottom,2)

    for i=1:nTop  

        if mode == "dependent"
            output[i,:,:] = generateDependentDirichletProcess(nBottom,alpha,probability)
        end

        if mode == "finite"
            output[i,:,:] = generateFiniteDirichletProcess(nBottom,alpha,probability)
        end 

        if mode == "stick"
            output[i,:,:] = generateStickBreaking(nBottom,alpha,probability)
        end


    end

    return output 



end
