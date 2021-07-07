using Images, ImageMorphology
using Statistics, StatsBase, Distances
using LinearAlgebra
using ProgressMeter
using OptimalTransport
using Tulip
using JuMP, Clp

# Lexicographical ordering with red, green and blue.
function h_rgb(img)
    return channelview(img)[1,:,:]+channelview(img)[2,:,:]/256+channelview(img)[3,:,:]/(256*256)
end

function h_red(img)
    return channelview(img)[1,:,:]
end

function h_mean(img)
    return reshape(mean(channelview(img),dims=1),size(img))
end

function marginalMM(img, mm_operator)
    img_m = channelview(copy(img))
    for i=1:3
        img_m[i,:,:] = mm_operator(img_m[i,:,:])
    end
    return colorview(RGB,img_m)
end

function hMM(img, h_function, mm_operator)
    n1,n2 = size(img)
    img_h = h_function(img)
    b = sortperm(img_h[:])
    im_latt = Array{Int64,1}(undef, n1*n2)
    im_latt[b] = 1:n1*n2
    im_latt = reshape(im_latt,n1,n2)
    img_mm = mm_operator(im_latt)
    return reshape(img[b[img_mm[:]]],n1,n2)
end

function img_hist(I)
    A = countmap(I)
    return float.(channelview(collect(keys(A)))), collect(values(A))
end
                       
function irregularity_distances(I,J,method="JuMP",epsilon=1.e-3,p=1)
    if p==2
        cost = sum((channelview(float.(I))-channelview(float.(J))).^2)
        metric = SqEuclidean()
    else
        cost = sum(sqrt.(sum((channelview(float.(I))-channelview(float.(J))).^2,dims=1)))
        metric = Euclidean()        
    end
    min_cost = cost
    if cost>0
        DI, HI = img_hist(I)
        DJ, HJ = img_hist(J)
        if (size(DI,2)>1)&&(size(DJ,2)>1)
            ND = prod(size(I))
            HI, HJ = HI/ND, HJ/ND
            C = pairwise(metric, DI, DJ, dims=2)
            if method=="sinkhorn_stabilized"
                # Solve the regularized transportation problem.
                plan = sinkhorn_stabilized(HI,HJ,C,epsilon)
                min_cost = ND*sum(C.*plan)
            elseif method=="sinkhorn_stabilized_epsscaling"
                plan = sinkhorn_stabilized_epsscaling(HI,HJ,C,epsilon)
                min_cost = ND*sum(C.*plan)
            elseif method=="sinkhorn"
                min_cost = ND*sinkhorn2(HI,HJ,C,epsilon,maxiter=5_000)
            elseif method=="emd"
                min_cost = ND*emd2(HI,HJ,convert(Array{Float64},C),Tulip.Optimizer())
            else
                # Solve the transportation problem using JuMP + Clp.
                NI, NJ = size(C)
                model = Model(Clp.Optimizer)
                @variable(model, X[1:NI,1:NJ]>=0)
                @objective(model, Min, sum(C.*X)) 
                @constraint(model, sum(X,dims=1).== reshape(HJ,(1,NJ)))
                @constraint(model, sum(X,dims=2).== reshape(HI,(NI,1)))
                set_silent(model)
                optimize!(model)
                min_cost = ND*objective_value(model)
            end
        end
    else 
        return [0,0]         
    end
    return [min_cost,cost]
end

function irregularity_global(I,J,method="JuMP",epsilon=1.e-3,p=1)
    if method=="JuMP"
        println("Global irregularity index computed analytically (with p = "*string(p)*").")
    else
        println("Global irregularity index computed using "*method*" (with epsilon = "*string(epsilon)*" and p = "*string(p)*")")
    end
    distances = irregularity_distances(I,J,method,epsilon,p)                    
    if distances[2]>0
        if p==2
            return 1-sqrt(distances[1]/distances[2])
        else
            return 1-distances[1]/distances[2]
        end 
    else 
        return 0         
    end
end                
                
function irregularity_local(I,J,Wsize = 16,method="sinkhorn_stabilized",epsilon=1.e-3,p=1)
    if method=="JuMP"
        println("Local irregularity index computed analytically (with p = "*string(p)*").")
    else
        println("Local irregularity index computed using "*method*" (with epsilon = "*string(epsilon)*" and p = "*string(p)*")")
    end
    ND = prod(size(I))
    M, N = div.(size(I),Wsize)
    println("Dimension of the local windows:", Wsize,"x",Wsize)
    distances = zeros(2)
                                        
    # Compute firts the last values!
    indI = 1+Wsize*(M-1):size(I)[1]
    for j=1:N-1
        indJ = 1+Wsize*(j-1):Wsize*j
        distances += irregularity_distances(I[indI,indJ],J[indI,indJ],method,epsilon,p)
    end
    indJ = 1+Wsize*(N-1):size(I)[2]
    distances += irregularity_distances(I[indI,indJ],J[indI,indJ],method,epsilon,p)  
    # Let's now compute the other values!
    @showprogress for i=1:M-1
        indI = 1+Wsize*(i-1):Wsize*i
        for j=1:N-1
            indJ = 1+Wsize*(j-1):Wsize*j
            distances += irregularity_distances(I[indI,indJ],J[indI,indJ],method,epsilon,p)
        end
        indJ = 1+Wsize*(N-1):size(I)[2]
        distances += irregularity_distances(I[indI,indJ],J[indI,indJ],method,epsilon,p)
    end
    if distances[2]>0
        if p==2
            return 1-sqrt(distances[1]/distances[2])
        else
            return 1-distances[1]/distances[2]
        end
    else
        return 0
    end
end

const irregularity_index = irregularity_local 
const global_irregularity = irregularity_global
const local_irregularity = irregularity_local
                                                
#
# Under development!
#
function irregularity_marginal(I,J,method="sinkhorn_stabilized",epsilon=1.e-3,p=2,metric=SqEuclidean())
    Ic = channelview(I)
    Jc = channelview(J)
    dists = zeros(2)
    for slice = 1:3
        cost = metric(Ic[slice,:,:],Jc[slice,:,:])
        min_cost = 0
        if cost>0
            ND = prod(size(I))
            DI, HI = img_hist(Ic[slice,:,:])
            DJ, HJ = img_hist(Jc[slice,:,:])
            if (length(DI)>1)&&(length(DJ)>1)
                HI, HJ = HI/ND, HJ/ND
                C = pairwise(metric, DI', DJ', dims=2)
                if method=="sinkhorn_stabilized"
                    # Solve the regularized transportation problem.
                    plan = sinkhorn_stabilized(HI,HJ,C,epsilon)
                    min_cost = ND*sum(C.*plan)
                elseif method=="emd"
                    min_cost = ND*emd2(HI,HJ,convert(Array{Float64},C),Tulip.Optimizer())
                else
                    min_cost = ND*sinkhorn2(HI,HJ,C,epsilon)
                end
            end
        end
        dists += [min_cost,cost]
    end
    if dists[2]>0
        if p==1
            return 1-dists[1]/dists[2]
        else
            return 1-sqrt(dists[1]/dists[2])
        end
    else
        return 0
    end
end
                                            
                                            
                                            
# function irregularity_index(I, J, Wstep = 4, verbose=true)
#     if verbose==true
#         println("Computing the irregularity using windows of size ",2*Wstep,"x",2*Wstep, ".")
#         n1,n2 = Int.(floor.(size(I)./Wstep))
#         # IrValues = zeros(n1-1,n2-1)
#         IrValue = 1.0
#         @showprogress for i=1:n1-1
#             for j=1:n2-1
#                 indI = 1+Wstep*(i-1):Wstep*(i+1)
#                 indJ = 1+Wstep*(j-1):Wstep*(j+1)
#                 # IrValues[i,j] = irregularity_small(I[indI,indJ],J[indI,indJ])
#                 IrValue = IrValue*(1-irregularity_global(I[indI,indJ],J[indI,indJ]))
#                 # print("\nIrregularity = ",IrValue)
#             end
#         end
#     else
#         n1,n2 = Int.(floor.(size(I)./Wstep))
#         # IrValues = zeros(n1-1,n2-1)
#         IrValue = 1.0
#         for i=1:n1-1
#             for j=1:n2-1
#                 indI = 1+Wstep*(i-1):Wstep*(i+1)
#                 indJ = 1+Wstep*(j-1):Wstep*(j+1)
#                 # IrValues[i,j] = irregularity_small(I[indI,indJ],J[indI,indJ])
#                 IrValue = IrValue*(1-irregularity_global(I[indI,indJ],J[indI,indJ]))
#                 # print("\nIrregularity = ",IrValue)
#             end
#         end
#     end
#     return 1-(IrValue)^(1/((n1-1)*(n2-1)))
# end
