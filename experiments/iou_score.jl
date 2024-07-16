using FileIO 
using JLD2
using ImageContrastAdjustment
using ImageView
using Images
using ImageFiltering
using PlotlyJS
using Plots

#Load in data
# file = "/Users/mward19/Documents/Segmentation/segment-bacteria/experiments/sample_data/run_6076.jld2"
# data = load_object(file)
# slice2 = data[100, :, :]

""" Subsample array by a factors which are powers of 2 greater than or equal to 1. """
function downsample(array, factors::AbstractVector)
    if length(factors) != ndims(array)
        throw(
            ErrorException("Each element of factors should correspond to an axis of the array.")
        )
    end
    for factor in factors
        if !ispow2(factor) || factor < 1
            throw(ErrorException("$factor is not a power of 2 greater than or equal to 2."))
        end
    end

    indices = [firstindex(i_range):factor:lastindex(i_range)
                for (i_range, factor) in zip(axes(array), factors)]
    return array[indices...]
end

segmentation = "/Users/mward19/Documents/Segmentation/segment-bacteria/experiments/sample_data/filled_seg_6074.jld2"
seg_data = downsample(load_object(segmentation),[8,8,8])


println(size(seg_data))
# datax = range(-150,150,length=300)
# datay = range(-480,480,length=960)
# dataz = range(-464,464,length=928)

datax = range(-19,19,length=38)
datay = range(-60,60,length=120)
dataz = range(-58,58,length=116)

X,Y,Z = mgrid(datax, datay, dataz)
# Prepare the value array (intensity or values at each point)
values = seg_data 

vol = volume(
    x=X[:],
    y=Y[:],
    z=Z[:],
    value=values[:],
    isomin=0.1,
    isomax=0.8,
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=17, # needs to be a large number for good volume rendering
)

plotly()
PlotlyJS.plot(vol)