import .Rays

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
    cc = Rays.closest_contour(𝐈, 𝐦, θ)

    display(contours)
    img = RGB.(𝐈.intensities)
    draw!(img, Cross(Point(reverse(𝐦)...), 30), RGB{Float64}(1,0,0))
    draw!(img, Cross(Point(reverse(cc)...), 30), RGB{Float64}(0,1,0)) |> display
end

function orientation_demo(𝐈, 𝐦, θ_target)
    # Make image smoother
    𝛉 = LinRange(0, 2π, 1000)
    orientations = [Rays.get_orientation(𝐈, 𝐦, θ) for θ in 𝛉]
    p = plot(
        𝛉 ./ π,
        orientations,
        xlabel = "π radians",
        ylabel = "Inner product with image gradient",
        legend=false
    )
    # Plot the point we are interested in
    scatter!(p, [θ_target/π], [Rays.get_orientation(𝐈, 𝐦, θ_target)])
    display(p)
end

function norm_demo(𝐈, 𝐦, θ_target)
    𝛉 = LinRange(0, 2π, 10000)
    norms = [Rays.get_norm(𝐈, 𝐦, θ) for θ in 𝛉]
    p = plot(
        𝛉 ./ π,
        norms,
        xlabel = "π radians",
        ylabel = "Gradient norm",
        legend=false
    )
    # Plot the point we are interested in
    scatter!(p, [θ_target/π], [Rays.get_norm(𝐈, 𝐦, θ_target)])
    display(p)
end

function dist_difference_demo(𝐈, 𝐦, θ_target)
    𝛉′ = LinRange(0, 2π, 1000) .+ θ_target
    dist_diffs = [Rays.get_dist_difference(𝐈, 𝐦, θ, θ′) for θ′ in 𝛉′]

    p = plot(
        (𝛉′ .- θ_target) ./ π ,
        dist_diffs,
        xlabel = "π radians from θ",
        ylabel = "Distance difference",
        legend=false
    )
    display(p)
end

𝐦 = [400, 400]
θ = π/2
𝐈 = Rays.Image(intensities, contours)
cc_demo(𝐈, 𝐦, θ)
#orientation_demo(𝐈, 𝐦, θ)
#norm_demo(𝐈, 𝐦, θ)
dist_difference_demo(𝐈, 𝐦, θ)