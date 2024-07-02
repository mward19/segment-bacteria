module FastRays
""" Implements optimizations in §3.3. """

using LinearAlgebra
using ImageFiltering
using DFT

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
            # TODO: implement
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

""" Calculates the gradient of 𝐈 at 𝐦 with Sobel kernels. Assumes 2D image. """
function imgrad(𝐈, 𝐦)
    Δ = 1 # Sobel kernel is 3×3, so we only pad a pixel by 1 in each direction
    sobel_x = [
        -1  0  1
        -2  0  2
        -1  0  1
    ]
    sobel_y = [
         1  2  1
         0  0  0
        -1 -2 -1
    ]
    # TODO: implement
end

end # module