module FastRays
""" Implements optimizations in §3.3. """

using LinearAlgebra
using ImageFiltering
using DSP

struct Image
    intensities::Array
    contours::Array
end

struct RayData
    𝐈::Image
    # The three precomputable Ray features
    dist::Array
    ori:: Array
    norm::Array
    calculated::Array{Bool}
end

RayData(𝐈) = RayData(
    𝐈, 
    Array{Float64}(undef, size(𝐈)...), 
    Array{Float64}(undef, size(𝐈)...), 
    Array{Float64}(undef, size(𝐈)...),
    fill(false, size(𝐈)...)
) 


struct ScanLine
    𝐈::Image
    𝐦::Vector # Initial position
    θ::Number # Angle
    𝐮::Vector # Step direction (unit vector)
    rd::RayData
end
ScanLine(𝐈, 𝐦, θ) = ScanLine(𝐈, 𝐦, θ, [cos(θ), sin(θ)], nothing)
Base.iterate(S::ScanLine, state=𝐦) = in_bounds(𝐈, state) ? (state + S.𝐮) : nothing


# this is all wrong
function precompute(𝐈::Image, rd::RayData, 𝐦, θ)
    𝐮 = [cos(θ), sin(θ)] # scan line direction
    S = ScanLine(𝐈, 𝐦, θ, 𝐮, rd)
    
    d = 0
    o = 1
    n = 0
    for location in S
        𝐦 = floor.(Int, location) 
        # If already calculated, don't calculate again
        if rd.calculated[𝐦...]
            continue
        end
        
        rd.dist[𝐦...] = d
        rd.ori[𝐦...]  = o
        rd.norm[𝐦...] = n

        if is_contour(𝐈, 𝐦)
            d = 0
            grad = imgrad(𝐈.intensities, 𝐦)
            n = norm(grad)
            o = 1/n * grad ⋅ 𝐮
        else
            d += 1
        end
    end
    
end


""" Checks that a given location 𝐦 is in the bounds of `image`. """
function in_bounds(𝐈::Image, 𝐦::Vector)
    return 𝐦 == loc_in_image(𝐈, 𝐦)
end

"""
If the given point is outside the image bounds, represent the out of bounds dimensions as ∞ or -∞.
"""
function loc_in_image(𝐈::Image, 𝐦::Vector)
    𝐦 = Float64.(𝐦)
    image_dims = size(𝐈.intensities)
    @assert length(image_dims) == length(𝐦)

    for i in eachindex(image_dims)
        if floor(𝐦[i]) < 1
            𝐦[i] = -Inf
        elseif floor(𝐦[i]) > image_dims[i]
            𝐦[i] = Inf
        end
    end
    return 𝐦
end

""" Checks if a given location 𝐦 is a contour in the image 𝐈. """
is_contour(𝐈::Image, 𝐦::Vector) = (0 != 𝐈.contours[floor.(Int, 𝐦)...])

""" Assists in providing mirror padding for an image, but by index. """
function mirror_index(intensities, 𝐦)
    new_𝐦 = 𝐦
    end_indices = size(intensities)
    for i in eachindex(𝐦)
        begin_index = 1
        end_index = end_indices[i]
        if 𝐦[i] < begin_index
            offset = begin_index - 𝐦[i]
            new_𝐦[i] = begin_index + offset
        elseif 𝐦[i] > end_index
            offset = end_index - 𝐦[i]
            new_𝐦[i] = end_index + offset
        end
    end
    return new_𝐦
end         

""" Get mirror padded (radius of Δ) neighborhood of image intensities around 𝐦. """
function neighborhood(intensities, 𝐦, Δ=1)
    neigh = Array{eltype(intensities)}(undef, fill(2Δ+1, ndims(intensities))...)
    offsets = Iterators.product(fill(-Δ:Δ, length(𝐦))...) # offsets ``
    to_index = fill(Δ+1, length(𝐦)) # offsets `offsets` to get 1-index values
    for offset in offsets
        offset = collect(offset)
        mirror_coord = mirror_index(intensities, 𝐦 + offset)
        neigh[(offset + to_index)...] = intensities[mirror_coord...]
    end

    return neigh
end


""" Calculates the gradient of intensities at 𝐦 with Sobel kernels. Assumes 2D image. """
function imgrad(intensities, 𝐦)
    sobel_x = [
        -1 -2 -1
         0  0  0
         1  2  1
    ] ./8
    sobel_y = [
        -1  0  1
        -2  0  2
        -1  0  1
    ] ./8

    Δ = 1 # Sobel kernel is 3×3, so we only pad a pixel by 1 in each direction
    kernel_center = fill(1 + Δ, length(𝐦))
    dim_ranges = fill(-Δ:Δ, length(𝐦))
    neigh = neighborhood(intensities, 𝐦, Δ)
    x_gradient = +((neigh .* sobel_x)...)
    y_gradient = +((neigh .* sobel_y)...)
    return [x_gradient, y_gradient]
end

end # module