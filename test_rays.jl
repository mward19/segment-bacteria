import .Rays2D

using Images, ImageEdgeDetection, ImageDraw
using ImageEdgeDetection: Percentile
using JLD2
using LocalFilters
using ImageFiltering
using LinearAlgebra
using Plots
using Infiltrator

intensities = Float64.(load_object("data/raw_slice_5.jld2"))

# Apply bilateral filters
filters = [
    img -> bilateralfilter(img, 1.0, 5), 
    img -> bilateralfilter(img, 40/256, 8),
    #img -> bilateralfilter(img, 20/256, 15)
]

# Apply filters and display results
for f in filters
    intensities = f(intensities)
end
intensities = Gray.(intensities)
# Done filtering.

canny(σ) = Canny(spatial_scale=σ, high=Percentile(98), low=Percentile(80))
contours = detect_edges(intensities, canny(5))

function cc_demo(𝐈, 𝐦, θ)
    cc = Rays2D.closest_contour(𝐈, 𝐦, θ)

    display(contours)
    img = RGB.(𝐈.intensities)
    draw!(img, Cross(Point(reverse(𝐦)...), 30), RGB{Float64}(1,0,0))
    draw!(img, Cross(Point(reverse(cc)...), 30), RGB{Float64}(0,1,0)) |> display
end

function orientation_demo(𝐈, 𝐦, θ_target)
    # Make image smoother
    intens = imfilter(𝐈.intensities, Kernel.gaussian(16))
    new_𝐈 = Rays2D.Image(intens, 𝐈.contours)
    𝚯 = LinRange(0, 2π, 10000)
    orientations = [Rays2D.orientation(new_𝐈, 𝐦, θ) for θ in 𝚯]
    p = plot(
        𝚯 ./ π,
        orientations,
        xlabel = "π radians",
        ylabel = "Inner product with image gradient",
        legend=false
    )
    # Plot the point we are interested in
    scatter!(p, [θ_target/π], [Rays2D.orientation(new_𝐈, 𝐦, θ_target)])
    display(p)
end

function norm_demo(𝐈, 𝐦, θ_target)
    𝚯 = LinRange(0, 2π, 10000)
    norms = [Rays2D.norm(𝐈, 𝐦, θ) for θ in 𝚯]
    p = plot(
        𝚯 ./ π,
        norms,
        xlabel = "π radians",
        ylabel = "Gradient norm",
        legend=false
    )
    # Plot the point we are interested in
    scatter!(p, [θ_target/π], [Rays2D.norm(𝐈, 𝐦, θ_target)])
    display(p)
end

# Testing orientation feature
𝐈 = Rays2D.Image(intensities, contours)
𝐦 = [445, 450]
θ = .15 * π
cc_demo(𝐈, 𝐦, θ)
#plot(intensities) |> display
orientation_demo(𝐈, 𝐦, θ)
norm_demo(𝐈, 𝐦, θ)

Rays2D.dist_difference(𝐈, 𝐦, 2π/3, π)