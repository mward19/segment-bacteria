# https://www.tugraz.at/fileadmin/user_upload/Institute/ICG/Images/team_lepetit/publications/smith_iccv09.pdf

module Rays2D

struct Image
    intensities::Array
    contours::Array
end

""" Checks that a given location 𝐦 is in the bounds of `image`. """
function in_bounds(𝐈::Image, 𝐦::Tuple)
    return 𝐦 == image_loc(𝐈, 𝐦)
end

"""
If the given point is outside the image bounds, represent the out of bounds dimensions as ∞ or -∞.
"""
function loc_in_image(𝐈::Image, 𝐦::Tuple)
    𝐦_new = 𝐦
    image_dims = size(𝐈.intensities)
    for (i, Nᵢ, 𝐦ᵢ) in enumerate(zip(image_dims, 𝐦))
        if floor(𝐦ᵢ) < 1
            𝐦_new[i] = -Inf
        elseif floor(𝐦ᵢ) > Nᵢ
            𝐦_new[i] = Inf
        end
    end
    return 𝐦_new
end

""" Checks if a given location 𝐦 is a contour in the image 𝐈. """
is_contour(𝐈::Image, 𝐦::Tuple) = (0 != 𝐈.contours[floor.(𝐦)...])

""" Closest contour point 𝐜. """
function closest_contour(𝐈::Image, 𝐦::Tuple, θ::Number)
    step = (cos(θ), sin(θ))
    step ./= norm(step) # Normalize step

    while in_bounds(𝐈, 𝐦)
        if is_contour(𝐈, 𝐦)
            return floor.(𝐦)
        end
        𝐦 .+= step
    end

    # If there is no contour before the edge of the image,
    # use `loc_in_image` to place ∞ values.
    return floor.(loc_in_image(𝐦))
end

end # module