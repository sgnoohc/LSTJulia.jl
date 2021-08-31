module LSTJulia

using CSV
using DataFrames
using Arrow
using Query
using PlotlyJSWrapper
using FHist
using ProgressMeter

export studyfishbone, writearrow, make_fb_mult_plots, make_angle_plots

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
    smd = sharemd(v1, v2)
    if smd == 1
        dx = v1.hit5x-v2.hit5x
        dy = v1.hit5y-v2.hit5y
        dxy = sqrt(dx^2+dy^2)
    else
        dx = v1.hit4x-v2.hit4x
        dy = v1.hit4y-v2.hit4y
        dxy = sqrt(dx^2+dy^2)
    end
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
function getangledata(evt_)
    segs = Dict{Int, Vector{Segment}}()
    nsegs = Int[]

    # println(length(evt_))

    evt = filter(x->istruth(x), evt_)

    # println(evt_ .|> x->x.simidx1)

    # Get all the possible simidxs
    simidxs = evt .|> x->x.simidx1

    # Loop over possible simidxs
    for simidx in Set(simidxs)

        # Aggregate the segments
        sgs = filter(x->x.simidx1 == simidx, evt)
        segs[simidx] = sgs
        push!(nsegs, length(sgs))

    end

    angledata = SegmentPairAngle[]
    for (simidx, sgs) in segs
        for i in 1:length(sgs)
            for j in i+1:length(sgs)
                smd = sharemd(sgs[i], sgs[j])
                moddiff = moduledifftype(sgs[i], sgs[j])
                Δϕ = deltaphi(sgs[i], sgs[j])
                Δθ = deltatheta(sgs[i], sgs[j])
                Δxy = deltaxy(sgs[i], sgs[j])
                if smd >= 1
                    push!(angledata, SegmentPairAngle(smd, moddiff, Δϕ, Δθ, Δxy))
                elseif smd < 0
                    println("ERROR: these are same md! Why are you comparing same segments?")
                end
            end
        end
    end
    return (angledata, nsegs)
end

# _________________________________________________________________________________________________________________________________
function studyfishbone(fn)

    # Full file name
    fulln = "data/segment/$fn.arrow"

    @info "Opening up file $fulln"

    # Open up the segments data
    table = Arrow.Table("$fulln")

    # All the segment pairs to be plotted
    sgpairs = Vector{SegmentPairAngle}()
    sgpairs_fb = Vector{SegmentPairAngle}()
    nsims = Vector{Int}()
    ndups = Vector{Int}()
    nsims_fb = Vector{Int}()
    ndups_fb = Vector{Int}()

    @info "Looping..."

    # Loop over the segment data
    @showprogress for evt in table.event
        # Regular
        (a, n) = getangledata(evt)
        append!(sgpairs, a)
        nsim = length(n)
        push!(nsims, nsim)
        ndup = length(filter(x->x>1, n))
        push!(ndups, ndup)
        # Run fishbone
        evt_fb = runfishbone(evt)
        # Fishboned
        (a, n) = getangledata(evt_fb)
        append!(sgpairs_fb, a)
        nsim_fb = length(n)
        push!(nsims_fb, nsim_fb)
        ndup_fb = length(filter(x->x>1, n))
        push!(ndups_fb, ndup_fb)
    end

    df = DataFrame(
                   nsims = nsims,
                   nsims_fb = nsims_fb,
                   ndups = ndups,
                   ndups_fb = ndups_fb,
                  )

    @time Arrow.write("data/angle/$fn.arrow", sgpairs)
    @info "Wrote data/angle/$fn.arrow"
    @time Arrow.write("data/anglefb/$fn.arrow", sgpairs_fb)
    @info "Wrote data/anglefb/$fn.arrow"
    @time Arrow.write("data/fb/$fn.arrow", df)
    @info "Wrote data/fb/$fn.arrow"

end

# # _________________________________________________________________________________________________________________________________
# function getangledata(fn)
#     table = Arrow.Table("data/angle/$fn.arrow")
#     df = DataFrame(table)
#     as = []
#     for item in 1:4
#         angledata = df |> @filter(_.moddifftype == item) |> DataFrame
#         push!(as, angledata)
#     end
#     return as
# end

# _________________________________________________________________________________________________________________________________
function make_fb_mult_plots(fn)

    table = Arrow.Table("data/fb/$fn.arrow")
    df = DataFrame(table)

    bins = if occursin("PU200", fn)
        0:20:500
    else
        0:1:20
    end

    nsims = Hist1D(df.nsims, bins)
    nsims_fb = Hist1D(df.nsims_fb, bins)

    bins = if occursin("PU200", fn)
        0:1:100
    else
        0:1:20
    end

    ndups = Hist1D(df.ndups, bins)
    ndups_fb = Hist1D(df.ndups_fb, bins)

    p = plot_stack(
         backgrounds=[nsims],
         signals=[nsims_fb],
         outputname=string("plots/",fn,"_nsims.{html,pdf}"),
         backgroundlabels=["before fishbone"],
         signallabels=["after fishbone"],
         xaxistitle="N<sub>sims</sub>",
         showtotallegend=false,
         showbeaminfo=false,
         backgroundcolors=[6004],
         cmsextralabeltext="Simulation",
        );

    p = plot_stack(
         backgrounds=[ndups],
         signals=[ndups_fb],
         outputname=string("plots/",fn,"_ndups.{html,pdf}"),
         backgroundlabels=["before fishbone"],
         signallabels=["after fishbone"],
         xaxistitle="N<sub>dups</sub>",
         showtotallegend=false,
         showbeaminfo=false,
         backgroundcolors=[6004],
         cmsextralabeltext="Simulation",
        );

    return

end

# _________________________________________________________________________________________________________________________________
function make_angle_plots(fn)
    table = Arrow.Table("data/angle/$fn.arrow")
    df = DataFrame(table)
    as = []
    for item in 1:4
        angledata = df |> @filter(_.moddifftype == item) |> DataFrame
        push!(as, angledata)
    end

    ylog = occursin("PU200", fn)
    yrange = [-1, 10]
    ylog = false
    yrange = [0, 500]

    hs = as .|> x->Hist1D(abs.(x.Δϕ), 0:0.0005:0.02)

    p = plot_stack(
                   backgrounds=[hs[4]],
                   signals=[hs[1], hs[2], hs[3]],
                   outputname=string("plots/",fn,"_deltaphi.{html,pdf}"),
                   backgroundlabels=["moddiff4"],
                   signallabels=["moddiff1", "moddiff2", "moddiff3"],
                   xaxistitle="Δϕ [rad]",
                   showtotallegend=false,
                   showbeaminfo=false,
                   backgroundcolors=[6004],
                   cmsextralabeltext="Simulation",
                   ylog=ylog,
                   yrange=yrange,
                  );

    hs = as .|> x->Hist1D(abs.(x.Δθ), 0:0.0005:0.02)

    p = plot_stack(
                   backgrounds=[hs[4]],
                   signals=[hs[1], hs[2], hs[3]],
                   outputname=string("plots/",fn,"_deltarz.{html,pdf}"),
                   backgroundlabels=["moddiff4"],
                   signallabels=["moddiff1", "moddiff2", "moddiff3"],
                   xaxistitle="Δrz [rad]",
                   showtotallegend=false,
                   showbeaminfo=false,
                   backgroundcolors=[6004],
                   cmsextralabeltext="Simulation",
                   ylog=ylog,
                   yrange=yrange,
                  );

    hs = as .|> x->Hist1D(abs.(x.Δxy), 0:0.05:4)

    p = plot_stack(
                   backgrounds=[hs[4]],
                   signals=[hs[1], hs[2], hs[3]],
                   outputname=string("plots/",fn,"_deltaxy.{html,pdf}"),
                   backgroundlabels=["moddiff4"],
                   signallabels=["moddiff1", "moddiff2", "moddiff3"],
                   xaxistitle="Δxy [cm]",
                   showtotallegend=false,
                   showbeaminfo=false,
                   backgroundcolors=[6004],
                   cmsextralabeltext="Simulation",
                   ylog=ylog,
                   yrange=yrange,
                  );


        # h = Hist1D(abs.(angledata.Δϕ), 0:0.0005:0.02)
        # p = plot_stack(
        #                backgrounds=[h],
        #                signals=[nsims_fb],
        #                outputname="nsims.{html,pdf}",
        #                backgroundlabels=["before fishbone"],
        #                signallabels=["after fishbone"],
        #                xaxistitle="N<sub>sims</sub>",
        #                showtotallegend=false,
        #                showbeaminfo=false,
        #               );



        # p = plot(trace, Layout(title="Δϕ distribution diffmode=$item", xaxis_title="Δϕ [rad]", yaxis_title="Events", bargap=0, template="simple_white"))
        # savefig(p, string("plots/moddifftypewide",item,"_",jobidx,"_dphi.pdf"))
        # savefig(p, string("plots/moddifftypewide",item,"_",jobidx,"_dphi.png"))

        # trace = PlotlyJSWrapper.build_hist1dtrace(PlotlyJSWrapper.make_fhist1d(abs.(angledata.Δθ), 0:0.0005:0.02), witherror=true)
        # # trace.fields[:marker][:color] = "blue"
        # # trace.fields[:opacity] = 0.75
        # p = plot(trace, Layout(title="Δrz distribution diffmode=$item", xaxis_title="Δrz [rad]", yaxis_title="Events", bargap=0, template="simple_white"))
        # savefig(p, string("plots/moddifftypewide",item,"_",jobidx,"_drz.pdf"))
        # savefig(p, string("plots/moddifftypewide",item,"_",jobidx,"_drz.png"))
    # end
end

end
