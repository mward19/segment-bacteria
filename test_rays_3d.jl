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

function display3d(orig_data, points=[]; rad=3)
    data = RGB.(orig_data)
    colors = [RGB(1,0,0), RGB(0,1,0), RGB(0, 0, 1)]
    Δ = rad
    for slice_index in axes(data, 1)
        for (p_index, p) in enumerate(points)
            if p[1] == slice_index
                data[p[1]-Δ:p[1]+Δ, p[2]-Δ:p[2]+Δ, p[3]-Δ:p[3]+Δ] .= colors[p_index]
                #data[p[1]-Δ:p[1]+Δ, p[2], p[3]] .= colors[p_index]
                #data[p[1], p[2]-Δ:p[2]+Δ, p[3]] .= colors[p_index]
                #data[p[1], p[2], p[3]-Δ:p[3]+Δ] .= colors[p_index]
            end
        end
        display(data[slice_index, :, :])
    end
end

original = rescale(load_object("data/mba2011-04-12-13.jld2"))
filters = [
    img -> bilateralfilter(img, 1.0, 5), 
    img -> bilateralfilter(img, 40/256, 8),
    #img -> bilateralfilter(img, 20/256, 15)
]

# Apply filters and display results
intensities = original
for f in filters
    intensities = f(intensities)
end
intensities = Float64.(intensities)
#display3d(intensities)


image = Rays.Image(intensities)