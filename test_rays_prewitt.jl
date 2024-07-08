import .Rays

using Images, ImageEdgeDetection, ImageDraw
using ImageEdgeDetection: Percentile
using JLD2
using LocalFilters
using ImageFiltering
using LinearAlgebra
using Plots
using DSP
using ArrayPadding
using StatsBase
using Infiltrator

function rescale(array)
    global_min = min(array...)
    global_max = max(array...)
    return (array .- global_min) ./ (global_max - global_min)
end

original = Float64.(load_object("data/mba2011-04-12-13.jld2")[100, :, :])
intensities = rescale(original)
# Apply bilateral filters
filters = [
    img -> bilateralfilter(img, 1.0, 5), 
    img -> bilateralfilter(img, 40/256, 8),
    img -> bilateralfilter(img, 20/256, 15)
]

# Apply filters and display results
for f in filters
    intensities = f(intensities)
end
intensities = Float64.(intensities)
# Done filtering.

# Calculate gradient image
prewitt = [
    [
        [1., 1., 1.], [-1., 0., 1.]
    ],
    [
        [-1., 0., 1.], [1., 1., 1.]
    ]
]
Δ = 1 # By how much the gradient image is padded by the convolution
# TODO: make more efficient with PaddedViews.jl
intensities_padded = pad(intensities, :mirror, Δ)
grad_1 = conv(prewitt[1][1], prewitt[1][2], intensities_padded)[begin+Δ:end-Δ,begin+Δ:end-Δ]
grad_2 = conv(prewitt[2][1], prewitt[2][2], intensities_padded)[begin+Δ:end-Δ,begin+Δ:end-Δ]
m, n = size(intensities_padded)
gradient = [[grad_1[i,j], grad_2[i,j]] for i in 1+Δ:m-Δ, j in 1+Δ:n-Δ]
gradient_norm = norm.(gradient)

# Calculate contours.
# TODO: using percentiles will cause problems in images with lots of edge noise (carbon holes)
percentile_threshold = 90
threshold = percentile(vec(gradient_norm), percentile_threshold)
contours = BitArray([gradient_norm[i...] >= threshold 
            for i in Iterators.product(axes(gradient_norm)...)])

image = Rays.Image(intensities, contours, gradient, gradient_norm)
𝐦 = [63, 300]
θ = - π/4 + 0.01
cc = Rays.closest_contour(image, 𝐦, θ)

show_pos(image, vec) = scatter!(plot(Gray.(image.contours)), [vec[2]], [vec[1]])
show_pos(image, 𝐦) |> display
show_pos(image, cc) |> display

orient = Rays.get_orientation(image, 𝐦, θ)
# Working!!