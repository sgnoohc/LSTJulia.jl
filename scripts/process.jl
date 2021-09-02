Revise.revise_file_now(Revise.pkgdatas[Revise.PkgId(LSTJulia)], "src/LSTJulia.jl")

using Revise
using LSTJulia
using PlotlyJSWrapper
using FHist

#----------------------------------------------

idxs = ["12", "23", "34", "45", "56"]

#----------------------------------------------

function make_plots_v1(hs, sample, idx, xaxis, name)
    h1 = hs[1] |> rebin(27)
    h2 = hs[2] |> rebin(27)
    h3 = hs[3] |> rebin(27)
    p = plot_stack(
                   backgrounds=[h3],
                   signals=[h1, h2],
                   outputname=string("plots/",sample,"_segments",idx,"_$name.{html,png,pdf}"),
                   backgroundlabels=["both"],
                   signallabels=["η stagger", "ϕ stagger"],
                   xaxistitle="$xaxis [rad]",
                   showtotallegend=false,
                   showbeaminfo=false,
                   backgroundcolors=[6004],
                   cmsextralabeltext="Simulation",
                  );
end

function make_mult_plots(hs, sample, idx)
    p = plot_stack(
         backgrounds=[hs[1]],
         signals=[hs[2]],
         outputname=string("plots/",sample,"_segments",idx,"_nsims.{html,png,pdf}"),
         backgroundlabels=["before fishbone"],
         signallabels=["after fishbone"],
         xaxistitle="N<sub>sims</sub>",
         showtotallegend=false,
         showbeaminfo=false,
         backgroundcolors=[6004],
         cmsextralabeltext="Simulation",
        );
    p = plot_stack(
         backgrounds=[hs[3]],
         signals=[hs[4]],
         outputname=string("plots/",sample,"_segments",idx,"_ndups.{html,png,pdf}"),
         backgroundlabels=["before fishbone"],
         signallabels=["after fishbone"],
         xaxistitle="N<sub>dups</sub>",
         showtotallegend=false,
         showbeaminfo=false,
         backgroundcolors=[6004],
         cmsextralabeltext="Simulation",
        );
end

function make_plots_v2(hists, sample, moddiff, var, name)
    myhs = []
    labels = []
    for idx in idxs
        h = (hists["$var$idx"] .|> rebin(27))[moddiff]
        push!(myhs, h)
        if moddiff == 1
            push!(labels, "η stagger $idx")
        elseif moddiff == 2
            push!(labels, "ϕ stagger $idx")
        elseif moddiff == 3
            push!(labels, "both stagger $idx")
        end
    end
    p = plot_stack(
                   backgrounds=[myhs[end]],
                   data=myhs,
                   outputname=string("plots/",sample,"_segments_",moddiff,"_$name","_log.{html,png,pdf}"),
                   backgroundlabels=[labels[end]],
                   datalabels=labels,
                   xaxistitle="$var [rad]",
                   showtotallegend=false,
                   showbeaminfo=false,
                   backgroundcolors=[6004],
                   cmsextralabeltext="Simulation",
                   signalcolors=repeat([8001, 8002, 8003, 8004, 6005],10),
                   ylog=true,
                   yrange=[0,5],
                  );
    p = plot_stack(
                   backgrounds=[myhs[end]],
                   signals=myhs,
                   outputname=string("plots/",sample,"_segments_",moddiff,"_$name","_lin.{html,png,pdf}"),
                   backgroundlabels=[labels[end]],
                   signallabels=labels,
                   xaxistitle="$var [rad]",
                   showtotallegend=false,
                   showbeaminfo=false,
                   backgroundcolors=[6004],
                   cmsextralabeltext="Simulation",
                   signalcolors=repeat([8001, 8002, 8003, 8004, 6005],10),
                  );
end

#---------------------------------------------

fns = []
for idx in idxs
    fn = "muon_gun_segments$idx"
    push!(fns, fn)
    fn = "muonhighe_gun_segments$idx"
    push!(fns, fn)
end

#----------------------------------------------

Threads.@threads for fn in fns
    @time writearrow(fn);
end

#----------------------------------------------

Threads.@threads for fn in fns
    @time studyfishbone(fn);
end

#----------------------------------------------

sample = "muon_gun"
hists = make_angle_hists(sample)
multhists = make_mult_hists(sample)

#----------------------------------------------

for idx in idxs
    make_plots_v1(hists["Δϕ$idx"], sample, idx, "Δϕ", "deltaphi")
    make_plots_v1(hists["Δθ$idx"], sample, idx, "Δθ", "deltarz")
    make_mult_plots(multhists["mult$idx"], sample, idx)
end
make_plots_v2(hists, sample, 1, "Δϕ", "deltaphi");
make_plots_v2(hists, sample, 2, "Δϕ", "deltaphi");
make_plots_v2(hists, sample, 3, "Δϕ", "deltaphi");
make_plots_v2(hists, sample, 1, "Δθ", "deltarz");
make_plots_v2(hists, sample, 2, "Δθ", "deltarz");
make_plots_v2(hists, sample, 3, "Δθ", "deltarz");
