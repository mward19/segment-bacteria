module Rays

using LinearAlgebra
using ImageFiltering
using StatsBase
using Infiltrator

struct Image
    intensities::AbstractArray
    contours::AbstractArray{Bool}
    normalized_gradient::AbstractArray # of gradient vectors
    gradient_norm::AbstractArray
    grad_calculated::AbstractArray{Bool} # true if gradient has been calculated
    cc_memo::Dict # Memo of closest contours. Key is tuple of angle(s)
end

#""" Constructs an Image from intensities and contours. """
#Image(intensities::AbstractArray, contours::AbstractArray) = Image(
#    intensities,
#    contours,
#    Array{Vector{Float64}}(undef, size(intensities)...),
#    Array{Float64}(undef, size(intensities)...),
#    AbstractArray{Bool}(fill(false, size(intensities)...)),
#    Dict{Tuple, Vector{Float64}}()
#)

""" Constructs an Image from intensities only. """
function Image(intensities::AbstractArray)
    # Calculate gradients
    gradient_by_dim = imgradients(intensities, KernelFactors.ando3, "reflect")
    gradient = [
        [gradient[𝐦...] for gradient in gradient_by_dim]
        for 𝐦 in Iterators.product(axes(intensities)...)
    ]
    gradient_norm = norm.(gradient)
    normalized_gradient = gradient ./ gradient_norm
    grad_calculated = trues(size(gradient))
    
    # TODO:
    # Use gradients to get contours with thresholding
    percentile_threshold = 98
    threshold = percentile(vec(gradient_norm), percentile_threshold)
    contours = BitArray([gradient_norm[i...] >= threshold 
                            for i in Iterators.product(axes(gradient_norm)...)])
    return Image(
        intensities,
        contours,
        normalized_gradient,
        gradient_norm,
        grad_calculated,
        Dict{Tuple, Vector{Float64}}()
    )
end


""" Constructs an Image from intensities, contours, gradients, and gradient norms. """
function Image(
    intensities::AbstractArray, 
    contours::AbstractArray{Bool}, 
    gradient::AbstractArray{<:AbstractVector}, 
    gradient_norm::AbstractArray
)
    normalized_gradient = gradient ./ gradient_norm
    grad_calculated = trues(size(intensities))
    cc_memo = Dict{Tuple, Vector{Float64}}()
    @infiltrate
    return Image(
        intensities,
        contours,
        normalized_gradient,
        gradient_norm,
        grad_calculated,
        cc_memo
    )
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

function ray_vector(θ, γ=nothing)
    if isnothing(γ) # 2D case
        return [cos(θ), sin(θ)]
    else # 3D case
        return [sin(γ), cos(θ)*cos(γ), sin(θ)*cos(γ)]
    end
end

function calc_grad(𝐈::Image, 𝐦::Vector)
    grad = imgrad(𝐈.intensities, 𝐦)
    𝐈.normalized_gradient[𝐦...] = grad ./ norm(grad)
    𝐈.gradient_norm[𝐦...] = norm(grad)
    𝐈.grad_calculated[𝐦...] = true
end

function get_normalized_grad(𝐈::Image, 𝐦::Vector)
    if !𝐈.grad_calculated[𝐦...]
        calc_grad(𝐈, 𝐦)
    end
    return 𝐈.normalized_gradient[𝐦...]
end

function get_grad_norm(𝐈::Image, 𝐦::Vector)
    if !𝐈.grad_calculated[𝐦...]
        calc_grad(𝐈, 𝐦)
    end
    return 𝐈.gradient_norm[𝐦...]
end

""" Closest contour point 𝐜. θ, γ in radians. """
function closest_contour(𝐈::Image, 𝐦::Vector, θ, γ=nothing)
    # TODO: test that this works
    if haskey(𝐈.cc_memo, (θ, γ))
        return 𝐈.cc_memo[(θ, γ)]
    end
    # Otherwise find it
    step = ray_vector(θ, γ) # Already normalized

    while in_bounds(𝐈, 𝐦)
        if is_contour(𝐈, 𝐦)
            return floor.(Int, 𝐦)
        end
        𝐦 += step
    end
    # If there is no contour before the edge of the image,
    # use `loc_in_image` to place ∞ values.
    cc = loc_in_image(𝐈, 𝐦)
    # Save and return
    𝐈.cc_memo[(θ, γ)] = cc
    return cc
end

""" Distance feature. """
function get_distance(𝐈::Image, 𝐦::Vector, θ, γ=nothing)
    𝐜 = closest_contour(𝐈, 𝐦, θ, γ)
    if (Inf in 𝐜) || (-Inf in 𝐜)
        return Inf
    end
    return LinearAlgebra.norm(𝐜 - 𝐦)
end

""" Orientation feature. """ # TODO: seems to have some issues. demo
function get_orientation(𝐈::Image, 𝐦::Vector, θ, γ=nothing)
    𝐜 = closest_contour(𝐈, 𝐦, θ, γ)
    if (Inf in 𝐜) || (-Inf in 𝐜)
        return NaN
    end
    return get_normalized_grad(𝐈, 𝐜) ⋅ ray_vector(θ, γ)
end

""" Norm feature. """
function get_norm(𝐈::Image, 𝐦::Vector, θ, γ=nothing)
    𝐜 = closest_contour(𝐈, 𝐦, θ, γ)
    if (Inf in 𝐜) || (-Inf in 𝐜)
        return NaN
    end
    return get_grad_norm(𝐈, 𝐜)
end

""" Distance difference feature. """
function get_dist_difference(
        𝐈 ::Image, 
        𝐦 ::Vector, 
        θ , 
        θ′, 
        γ =nothing, 
        γ′=nothing
    )
    𝐜  = closest_contour(𝐈, 𝐦, θ , γ )
    𝐜′ = closest_contour(𝐈, 𝐦, θ′, γ′)
    return (norm(𝐜 - 𝐦) - norm(𝐜′ - 𝐦)) / norm(𝐜 - 𝐦)
end

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

# TODO: make gradient more efficient
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