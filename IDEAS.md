# Canonical orientation
In the paper, the canonical orientation is the direction of greatest variance. In our case, it should probably be max orthogonal variance in the xy plane, since there are so many holes otherwise.

# Features to find in superpixels
- Rays, except gradient magnitude? or maybe all Rays
- local histogram thing (dr. morse)
- avg intensity
- a bunch of matlab imageprops, basically

# Segmentation
Segment the interior (including the membrane) of each bacteria in the tomogram, including extracellular vescicles(?). Determine the human-perceivable "top" and "bottom" of the image of the bacteria, and stop segmenting there.

# Evaluating filters
Apply filter, apply edge detector, do IoU of edge over the edge of our segmentations.