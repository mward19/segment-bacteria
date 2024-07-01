# https://www.tugraz.at/fileadmin/user_upload/Institute/ICG/Images/team_lepetit/publications/smith_iccv09.pdf

module Rays2D

using Infiltrator
using LinearAlgebra

struct Image
    intensities::Array
    contours::Array
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
    step /= norm(step) # Normalize step

    @infiltrate
    while in_bounds(𝐈, 𝐦)
        @infiltrate 𝐦[1] == 621
        if is_contour(𝐈, 𝐦)
            return floor.(Int, 𝐦)
        end
        𝐦 += step
        @infiltrate 𝐦[1] >= 620
    end

    # If there is no contour before the edge of the image,
    # use `loc_in_image` to place ∞ values.
    return floor.(Int, loc_in_image(𝐈, 𝐦))
end

""" Distance feature. """
function distance(𝐈::Image, 𝐦::Vector, θ::Number)
    𝐜 = closest_contour(𝐈, 𝐦, θ)
    return norm(𝐜 - 𝐦)
end

end # module