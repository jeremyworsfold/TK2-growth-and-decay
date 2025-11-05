if isdefined(@__MODULE__, :LanguageServer)
    include("../src/TK2.jl")
    using .TK2
else
    using TK2
end

include("plotting/plotting.jl")
theme = include("plotting/theme.jl")
theme isa Attributes && set_theme!(theme)
CairoMakie.activate!()

SCRIPTNAME = "04-LNA"
SIMULATE = false
NUM_SIMS = 100
SAVE = true
SAVEFIG = true

#################################################
BETAS = [0.0, 0.2]
base_config = load(SimParams, SCRIPTNAME)
configs = [@set base_config.params.beta = β for β in BETAS]
fnames = ["$SCRIPTNAME-sym", "$SCRIPTNAME-asym"]

f1 = fig_in_cm(9, 8.6, backgroundcolor=:transparent)
axs = [
    Axis(f1[1, 1], ylabel=L"P^\ast(s)", xticklabelsvisible=false, backgroundcolor=:transparent),
    Axis(f1[2, 1], ylabel=L"P^\ast(s)", xlabel=L"s", backgroundcolor=:transparent),
]
linkxaxes!(axs...)
rowgap!(f1.layout, Relative(0.02))

for (i, (ax, fname, conf)) in enumerate(zip(axs, fnames, configs))
    if SIMULATE
        data, conf = ensemble_custom_sim(conf, NUM_SIMS)
        if SAVE
            save_res_and_conf(fname, data, conf)
        end
    else
        data, conf = load_res_and_conf(fname)
    end

    #################################################
    @unpack N, kappa, beta = conf.params
    approx = from_params(PDMPStationarySol, kappa, beta)

    ϕ1, ϕ2 = TK2.bounds(approx)
    x = collect(range(ϕ1 - 3 / sqrt(N), ϕ2 + 4 / sqrt(N), 2000))

    Πp, Πm, ϕ = theoretical_pdfs(approx)

    ΠpLNA, ΠmLNA = theoretical_LNA_pdfs(approx, x, N)
    Π0LNA = ΠpLNA .+ ΠmLNA

    manual_hist!(ax, data["edges"], data["x1"]; color=COLORS[1], label=L"P_+^*(s)")
    manual_hist!(ax, data["edges"], data["x2"]; color=COLORS[2], label=L"P_-^*(s)")
    plot_theory!(ax, x, ΠpLNA, var=:x1, label=L"P_+^*(s)")
    plot_theory!(ax, x, ΠmLNA, var=:x2, label=L"P_-^*(s)")
    plot_theory!(ax, x, Π0LNA, var=:s, linestyle=:dot, label=L"P_0^*(s)")

    xlims!(ax, extrema(x))
    ylims!(ax, 0, maximum(Π0LNA) * 1.05)

    if i == 2
        axislegend(
            ax,
            padding=0,
            rowgap=2,
            patchsize=(10, 10),
            patchlabelgap=5,
            merge=true,
        )
    end
end

f1

if SAVEFIG
    save_fig_and_conf(SCRIPTNAME, f1, configs[2])
end
