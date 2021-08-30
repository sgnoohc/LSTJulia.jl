module LSTJulia

using CSV
using DataFrames
using Arrow

export studyfishbone, writearrow

# _________________________________________________________________________________________________________________________________
struct Segment
    hit0x::Float64
    hit0y::Float64
    hit0z::Float64
    hit1x::Float64
    hit1y::Float64
    hit1z::Float64
    hit2x::Float64
    hit2y::Float64
    hit2z::Float64
    hit3x::Float64
    hit3y::Float64
    hit3z::Float64
    hit4x::Float64
    hit4y::Float64
    hit4z::Float64
    hit5x::Float64
    hit5y::Float64
    hit5z::Float64
    idx0::Int64
    idx1::Int64
    idx2::Int64
    idx3::Int64
    idx4::Int64
    idx5::Int64
    simidx0::Int64
    simidx1::Int64
    simidx2::Int64
    simidx3::Int64
    simidx4::Int64
    simidx5::Int64
    isFlatBarrel0::Int64
    isFlatBarrel1::Int64
    isPS0::Int64
    isPS1::Int64
    detId0::Int64
    detId1::Int64
    layer0::Int64
    layer1::Int64
    moduleId0::Int64
    moduleId1::Int64
    rod0::Int64
    rod1::Int64
end
ArrowTypes.arrowname(::Type{Segment}) = :Segment
ArrowTypes.JuliaType(::Val{:Segment}) = Segment

# _________________________________________________________________________________________________________________________________
struct SegmentPairAngle
    mdsharetype::Int64
    moddifftype::Int64
    Δϕ::Float64
    Δθ::Float64
    Δxy::Float64
end
ArrowTypes.arrowname(::Type{SegmentPairAngle}) = :SegmentPairAngle
ArrowTypes.JuliaType(::Val{:SegmentPairAngle}) = SegmentPairAngle

# _________________________________________________________________________________________________________________________________
@inline function atan2(y, x)
    if x != 0
        return atan(y, x)
    end
    if y == 0
        return 0
    end
    if y > 0
        return π / 2
    end
    return -π / 2
end

# _________________________________________________________________________________________________________________________________
@inline function phimpipi(x::Float64)
    while x >= π
        x -= 2. * π;
    end
    while x < -π
        x += 2. * π;
    end
    return x
end

# _________________________________________________________________________________________________________________________________
@inline function theta(v1::Segment)
    dx = v1.hit5x-v1.hit4x
    dy = v1.hit5y-v1.hit4y
    dz = v1.hit5z-v1.hit4z
    drt = sqrt(dx^2 + dy^2)
    theta = phimpipi(π + atan2(-drt, -dz));
    return theta
end

# _________________________________________________________________________________________________________________________________
@inline function phi(v1::Segment)
    dx = v1.hit5x-v1.hit4x
    dy = v1.hit5y-v1.hit4y
    phi = phimpipi(π + atan2(-dy, -dx));
end

# _________________________________________________________________________________________________________________________________
@inline function deltatheta(v1::Segment, v2::Segment)
    dx1 = v1.hit5x-v1.hit4x
    dy1 = v1.hit5y-v1.hit4y
    dz1 = v1.hit5z-v1.hit4z
    drt1 = sqrt(dx1^2 + dy1^2)
    dx2 = v2.hit5x-v2.hit4x
    dy2 = v2.hit5y-v2.hit4y
    dz2 = v2.hit5z-v2.hit4z
    drt2 = sqrt(dx2^2 + dy2^2)
    arg = (dz2 * dz1 + drt2 * drt1) / sqrt(drt2 * drt2 + dz2 * dz2) / sqrt(drt1 * drt1 + dz1 * dz1);
    if arg >= 1
        return 0
    else
        return acos(arg)
    end
end

# _________________________________________________________________________________________________________________________________
@inline function deltaphi(v1::Segment, v2::Segment)
    return phimpipi(phi(v1) - phi(v2))
end

# _________________________________________________________________________________________________________________________________
@inline function deltaxy(v1::Segment, v2::Segment)
    dx = v1.hit5x-v2.hit5x
    dy = v1.hit5y-v2.hit5y
    dxy = sqrt(dx^2+dy^2)
end

# _________________________________________________________________________________________________________________________________
@inline function r3(v1::Segment)
    dx1 = v1.hit5x-v1.hit4x
    dy1 = v1.hit5y-v1.hit4y
    dz1 = v1.hit5z-v1.hit4z
    return sqrt(dx1^2 + dy1^2 + dz1^2)
end

# _________________________________________________________________________________________________________________________________
@inline function istruth(v1::Segment)
    v1.simidx0 != v1.simidx1 && return false
    v1.simidx0 != v1.simidx2 && return false
    v1.simidx0 != v1.simidx3 && return false
    return true
end

# _________________________________________________________________________________________________________________________________
function writearrow(fn)

    csv = CSV.File("data/csv/$fn.csv"; header=false)
    df = DataFrame(csv)
    gdf = groupby(df, ["Column1"])

    events = Vector{Segment}[]
    for i in 1:length(gdf)
        segs = Segment[]
        for row in eachrow(gdf[i])
            tdf = DataFrame(row)
            tdf = tdf[:,2:ncol(tdf)]
            v = [ tdf[1, i] for i in 1:ncol(tdf) ]
            push!(segs, Segment(v...))
        end
        push!(events, segs)
    end

    table = (event=events,)
    Arrow.write("data/segment/$fn.arrow", table)

end

# _________________________________________________________________________________________________________________________________
function sharemd(v1::Segment, v2::Segment)
    islowersame = v1.idx0 == v2.idx0 && v1.idx1 == v2.idx1
    isuppersame = v1.idx2 == v2.idx2 && v1.idx3 == v2.idx3
    if (!islowersame && isuppersame)
        return 2
    elseif (islowersame && !isuppersame)
        return 1
    elseif (!islowersame && !isuppersame)
        return 0
    elseif (islowersame && isuppersame)
        return -1
    end
end

# _________________________________________________________________________________________________________________________________
function moduledifftype(v1::Segment, v2::Segment)
    smd = sharemd(v1, v2)
    if smd == 2
        if (v1.moduleId0 != v2.moduleId0 && v1.rod0 == v2.rod0)
            return 1
        elseif (v1.rod0 != v2.rod0 && v1.moduleId0 == v2.moduleId0)
            return 2
        elseif (v1.moduleId0 != v2.moduleId0 && v1.rod0 != v2.rod0)
            return 3
        elseif (v1.moduleId0 == v2.moduleId0 && v1.rod0 == v2.rod0)
            return 4
        end
    elseif smd == 1
        if (v1.moduleId1 != v2.moduleId1 && v1.rod1 == v2.rod1)
            return 1
        elseif (v1.rod1 != v2.rod1 && v1.moduleId1 == v2.moduleId1)
            return 2
        elseif (v1.moduleId1 != v2.moduleId1 && v1.rod1 != v2.rod1)
            return 3
        elseif (v1.moduleId1 == v2.moduleId1 && v1.rod1 == v2.rod1)
            return 4
        end
    else
        return 0
    end
end

# _________________________________________________________________________________________________________________________________
function issegmentofinterest(v1::Segment)
    # Segment selection
    v1.layer0 != 1 && return false
    v1.layer1 != 2 && return false
    v1.isFlatBarrel0 != 1 && return false
    v1.isFlatBarrel1 != 1 && return false
    return true
end

# _________________________________________________________________________________________________________________________________
function add_segment!(segsdict::Dict{Int, Vector{Segment}}, nsegsdict::Dict{Int, Int}, v1::Segment)
    # Organize by simidx
    if haskey(segsdict, v1.simidx1)
        push!(segsdict[v1.simidx1], v1)
        nsegsdict[v1.simidx1] += 1
    else
        segsdict[v1.simidx1] = Segment[v1]
        nsegsdict[v1.simidx1] = 1
    end
end

# _________________________________________________________________________________________________________________________________
function organize_segments!(segments_dict::Dict{Int, Vector{Segment}}, nsegments_dict::Dict{Int, Int}, evt::Vector{Segment})
    # Loop over segments and organize by simidx
    for i in 1:length(evt)
        seg = evt[i]
        # !issegmentofinterest(seg) && continue
        !istruth(seg) && continue
        add_segment!(segments_dict, nsegments_dict, seg)
    end
end

# _________________________________________________________________________________________________________________________________
function get_sgpairangles(segments_dict)
    angledata = SegmentPairAngle[]
    for (simidx, segs) in segments_dict
        for i in 1:length(segs)
            for j in i+1:length(segs)
                smd = sharemd(segs[i], segs[j])
                moddiff = moduledifftype(segs[i], segs[j])
                Δϕ = deltaphi(segs[i], segs[j])
                Δθ = deltatheta(segs[i], segs[j])
                Δxy = deltaxy(segs[i], segs[j])
                if smd >= 1
                    push!(angledata, SegmentPairAngle(smd, moddiff, Δϕ, Δθ, Δxy))
                elseif smd < 0
                    println("ERROR: these are same md! Why are you comparing same segments?")
                end
            end
        end
    end
    return angledata
end

# _________________________________________________________________________________________________________________________________
function fishboneresult(v1::Segment, v2::Segment)
    moddiff = moduledifftype(v1, v2)
    moddiff == 0 || moddiff == 4 && return 0 # keep both
    Δϕ = abs(deltaphi(v1, v2))
    Δθ = abs(deltatheta(v1, v2))
    r3v1 = r3(v1)
    r3v2 = r3(v2)
    overlap = false

    # Boundary obtained with low pt muon gun sample (DEPRECATED)
    # if moddiff == 1 && Δϕ > 0.001 && Δϕ < 0.004 && Δθ < 0.01
    #     overlap = true
    # elseif moddiff == 2 && Δϕ > 0.005 && Δϕ < 0.016 && Δθ < 0.014
    #     overlap = true
    # elseif moddiff == 3 && Δϕ > 0.004 && Δϕ < 0.02  && Δθ < 0.014
    #     overlap = true
    # end

    # Boundary considering high muon gun as well
    if moddiff == 1 && Δϕ < 0.004 && Δθ < 0.014
        overlap = true
    elseif moddiff == 2 && Δϕ < 0.016 && Δθ < 0.014
        overlap = true
    elseif moddiff == 3 && Δϕ < 0.02  && Δθ < 0.014
        overlap = true
    end

    if overlap
        if r3v1 <= r3v2
            return 1
        else
            return 2
        end
    else
        return 0
    end
end

# _________________________________________________________________________________________________________________________________
function runfishbone(segs)
    evt = Vector{Segment}()
    for i in 1:length(segs)
        toremove = false
        for j in 1:length(segs)
            if fishboneresult(segs[i], segs[j]) == 2
                toremove = true
            end
        end
        !toremove && push!(evt, segs[i])
    end
    return evt
end

# _________________________________________________________________________________________________________________________________
function studyfishbone(fn)

    # Open up the segments data
    table = Arrow.Table("data/segment/$fn.arrow")

    # All the segment pairs to be plotted
    sgpairs = Vector{SegmentPairAngle}()
    sgpairs_fb = Vector{SegmentPairAngle}()

    # Loop over the segment data
    for evt in table.event

        df = DataFrame(evt)

        gdf = groupby(df, ["simidx1"])

        println(gdf)

        break

        # # information to process for the event
        # segs = Dict{Int, Vector{Segment}}()

        # # Loop over segments and organize by simidx
        # for i in 1:length(evt)

        #     # Get the segment
        #     seg = evt[i]

        #     # Check if it is a true segment
        #     !istruth(seg) && continue

        #     # Organize by simidx
        #     if haskey(segsdict, seg.simidx1)
        #         push!(segsdict[seg.simidx1], seg)
        #         nsegsdict[seg.simidx1] += 1
        #     else
        #         segsdict[seg.simidx1] = Segment[seg]
        #         nsegsdict[seg.simidx1] = 1
        #     end
        # end
    end

end

# _________________________________________________________________________________________________________________________________
function studyfishbone_v1(fn)

    # Open up the segments data
    table = Arrow.Table("data/segment/$fn.arrow")

    sgpairs = SegmentPairAngle[]
    sgpairs_fb = SegmentPairAngle[]

    # Loop over the events
    for evt in table.event
        # information to process for the event
        segments_dict = Dict{Int, Vector{Segment}}()
        nsegments_dict = Dict{Int, Int}()
        organize_segments!(segments_dict, nsegments_dict, evt)
        append!(sgpairs, get_sgpairangles(segments_dict))

        evt_fishboned = runfishbone(evt)
        segments_fishboned_dict = Dict{Int, Vector{Segment}}()
        nsegments_fishboned_dict = Dict{Int, Int}()
        organize_segments!(segments_fishboned_dict, nsegments_fishboned_dict, evt_fishboned)
        append!(segmentpairangledata_fishboned[Threads.threadid()], get_sgpairangles(segments_fishboned_dict))
        nsims = length(nsegments_dict)
        nsims_fishboned = length(nsegments_fishboned_dict)

        ndups = 0
        for (key, nseg) in nsegments_dict
            if nseg > 1
                ndups += 1
            end
        end
        ndups_fishboned = 0
        for (key, nseg) in nsegments_fishboned_dict
            if nseg > 1
                ndups_fishboned += 1
            end
        end
        # println("$nsims $nsims_fishboned $ndups $ndups_fishboned")
        push!(nsims_list, nsims)
        push!(nsims_fishboned_list, nsims_fishboned)
        push!(ndups_list, ndups)
        push!(ndups_fishboned_list, ndups_fishboned)
    end

    segmentpairangledata = Vector{SegmentPairAngle}[]
    for i in 1:Threads.nthreads()
        push!(segmentpairangledata, SegmentPairAngle[])
    end

    segmentpairangledata_fishboned = Vector{SegmentPairAngle}[]
    for i in 1:Threads.nthreads()
        push!(segmentpairangledata_fishboned, SegmentPairAngle[])
    end

    nsims_list = Int[]
    nsims_fishboned_list = Int[]
    ndups_list = Int[]
    ndups_fishboned_list = Int[]

    # Loop over the events
    @time for evt in table.event
        # information to process for the event
        segments_dict = Dict{Int, Vector{Segment}}()
        nsegments_dict = Dict{Int, Int}()
        organize_segments!(segments_dict, nsegments_dict, evt)
        append!(segmentpairangledata[Threads.threadid()], get_sgpairangles(segments_dict))

        evt_fishboned = runfishbone(evt)
        segments_fishboned_dict = Dict{Int, Vector{Segment}}()
        nsegments_fishboned_dict = Dict{Int, Int}()
        organize_segments!(segments_fishboned_dict, nsegments_fishboned_dict, evt_fishboned)
        append!(segmentpairangledata_fishboned[Threads.threadid()], get_sgpairangles(segments_fishboned_dict))
        nsims = length(nsegments_dict)
        nsims_fishboned = length(nsegments_fishboned_dict)

        ndups = 0
        for (key, nseg) in nsegments_dict
            if nseg > 1
                ndups += 1
            end
        end
        ndups_fishboned = 0
        for (key, nseg) in nsegments_fishboned_dict
            if nseg > 1
                ndups_fishboned += 1
            end
        end
        # println("$nsims $nsims_fishboned $ndups $ndups_fishboned")
        push!(nsims_list, nsims)
        push!(nsims_fishboned_list, nsims_fishboned)
        push!(ndups_list, ndups)
        push!(ndups_fishboned_list, ndups_fishboned)
    end

    df = DataFrame(
                   nsims = nsims_list,
                   nsims_fishboned = nsims_fishboned_list,
                   ndups = ndups_list,
                   ndups_fishboned = ndups_fishboned_list
                  )

    segmentpairangledata = vcat(segmentpairangledata...)
    segmentpairangledata_fishboned = vcat(segmentpairangledata_fishboned...)

    @time Arrow.write("data/angle/$fn.arrow")
    @time Arrow.write("data/anglefb/$fn.arrow")
    @time Arrow.write("data/fb/$fn.arrow")

end

end
