import .Rays2D

using Images, ImageEdgeDetection, ImageDraw
using ImageEdgeDetection: Percentile
using JLD2

intensities = load_object("data/raw_slice_5.jld2")
canny(σ) = Canny(spatial_scale=σ, high=Percentile(98), low=Percentile(80))
contours = detect_edges(intensities, canny(5))

𝐈 = Rays2D.Image(intensities, contours)
point = [480, 480]
cc = Rays2D.closest_contour(𝐈, point, π/6)

display(contours)
img = RGB.(𝐈.intensities)
draw!(img, Cross(Point(reverse(point)...), 30), RGB{Float64}(1,0,0))
draw!(img, Cross(Point(reverse(cc)...), 30), RGB{Float64}(0,1,0))

Rays2D.distance(𝐈, point, π/6)