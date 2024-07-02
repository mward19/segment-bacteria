# https://www.tugraz.at/fileadmin/user_upload/Institute/ICG/Images/team_lepetit/publications/smith_iccv09.pdf

module Rays2D

using Infiltrator
using LinearAlgebra
using ImageFiltering

struct Image
    intensities::Array
    contours::Array
    gradients::Array # of gradient vectors
end

""" Constructs an Image from just intensities and contours. """
function Image(intensities::Array, contours::Array)
    all_gradients = imgradients(Float64.(intensities), KernelFactors.sobel, "symmetric")
    gradients = Array{Vector}(undef, size(intensities)...)
    for 𝐦 in Iterators.product(axes(intensities)...)
        gradients[𝐦...] = [g[𝐦...] for g in all_gradients]
    end
    return Image(intensities, contours, gradients)
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

""" Closest contour point 𝐜. θ in radians."""
function closest_contour(𝐈::Image, 𝐦::Vector, θ::Number)
    step = [cos(θ), sin(θ)]
    step /= LinearAlgebra.norm(step) # Normalize step

    while in_bounds(𝐈, 𝐦)
        if is_contour(𝐈, 𝐦)
            return floor.(Int, 𝐦)
        end
        𝐦 += step
    end

    # If there is no contour before the edge of the image,
    # use `loc_in_image` to place ∞ values.
    return loc_in_image(𝐈, 𝐦)
end

""" Distance feature. """
function distance(𝐈::Image, 𝐦::Vector, θ::Number)
    𝐜 = closest_contour(𝐈, 𝐦, θ)
    if (Inf in 𝐜) || (-Inf in 𝐜)
        return Inf
    end
    return LinearAlgebra.norm(𝐜 - 𝐦)
end

""" Orientation feature. """ # TODO: seems to have some issues. demo
function orientation(𝐈::Image, 𝐦::Vector, θ::Number)
    𝐜 = closest_contour(𝐈, 𝐦, θ)
    if (Inf in 𝐜) || (-Inf in 𝐜)
        return NaN
    end
    return normalize(𝐈.gradients[𝐜...]) ⋅ [cos(θ), sin(θ)]
end

""" Norm feature. """
function norm(𝐈::Image, 𝐦::Vector, θ::Number)
    𝐜 = closest_contour(𝐈, 𝐦, θ)
    if (Inf in 𝐜) || (-Inf in 𝐜)
        return NaN
    end
    @infiltrate
    return LinearAlgebra.norm(𝐈.gradients[𝐜...])
end

""" Distance difference feature. """
function dist_difference(𝐈::Image, 𝐦::Vector, θ::Number, θ′::Number)
    𝐜  = closest_contour(𝐈, 𝐦, θ )
    𝐜′ = closest_contour(𝐈, 𝐦, θ′)
    lanorm = LinearAlgebra.norm
    return (lanorm(𝐜 - 𝐦) - lanorm(𝐜′ - 𝐦)) / lanorm(𝐜 - 𝐦)
end

end # module