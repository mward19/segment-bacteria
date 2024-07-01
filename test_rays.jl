import .Rays2D

using Images, ImageEdgeDetection
using ImageEdgeDetection: Percentile
using JLD2

intensities = load_object("data/raw_slice_5.jld2")
canny(σ) = Canny(spatial_scale=σ, high=Percentile(98), low=Percentile(80))
contours = detect_edges(intensities, canny(5))

Gray.(intensities) |> display
Gray.(contours) |> display