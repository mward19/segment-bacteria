module Tracers

struct Tracer
    resolution::Vector
    loc_begin::Vector
    θ::Number
    loc_end::Vector
end
function Tracer(resolution, loc_begin, θ)
    # TODO: calculate loc_end.
end

end # module