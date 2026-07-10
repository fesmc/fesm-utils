# make_figures.jl
# ---------------------------------------------------------------------------
# Regenerate all documentation figures for the fesm-utils Quarto site.
#
#   cd docs/figures
#   julia --project=. make_figures.jl
#
# The tsgen curves below re-implement the analytic (time-driven) forcing
# formulas from src/tsgen.f90 (function eval_open) purely for illustration —
# they are NOT called from the Fortran library. Keep them in sync if the
# analytic methods change. Outputs are committed PNGs so the site builds
# without Julia.
# ---------------------------------------------------------------------------

using CairoMakie
CairoMakie.activate!(px_per_unit = 2)

const OUT = @__DIR__

# fesm-utils / cosmo-blue palette
const BLUE   = colorant"#1f6feb"
const ORANGE = colorant"#e8710a"
const GREEN  = colorant"#2ca02c"
const GREY   = colorant"#57606a"
const RED     = colorant"#d1005c"

# A clean, consistent look for every figure
set_theme!(Theme(
    fontsize = 15,
    Axis = (
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = false,
        rightspinevisible = false,
        xlabelfont = :bold,
        ylabelfont = :bold,
    ),
))

# ---------------------------------------------------------------------------
# tsgen analytic forcing (mirror of eval_open in src/tsgen.f90)
# ---------------------------------------------------------------------------

f_start(fmin, fmax, dfsign) = dfsign > 0 ? fmin : fmax

function f_ramp_time(t; dt_init, dt_ramp, fmin, fmax, dfsign = 1.0)
    τ  = max(0.0, t - dt_init)
    fs = f_start(fmin, fmax, dfsign)
    rate = abs(fmax - fmin) / dt_ramp
    clamp(fs + dfsign * rate * τ, fmin, fmax)
end

function f_ramp_slope(t; dt_init, df_dt_max, fmin, fmax, dfsign = 1.0)
    τ  = max(0.0, t - dt_init)
    fs = f_start(fmin, fmax, dfsign)
    clamp(fs + dfsign * df_dt_max * τ, fmin, fmax)
end

function f_ramp_time_step(t; dt_init, dt_ramp, dt_conv, fmin, fmax, fconv, dfsign = 1.0)
    τ   = max(0.0, t - dt_init)
    fs  = f_start(fmin, fmax, dfsign)
    fext = dfsign > 0 ? fmax : fmin
    if τ <= dt_ramp
        fs + (fext - fs) * (τ / dt_ramp)
    elseif τ <= dt_ramp + dt_conv
        fext + (fconv - fext) * ((τ - dt_ramp) / dt_conv)
    else
        fconv
    end
end

function f_sin(t; dt_ramp, fmin, fmax)
    amp = 0.5 * (fmax - fmin)
    off = 0.5 * (fmax + fmin)
    amp * sin(2π * t / dt_ramp) + off
end

# ---------------------------------------------------------------------------
# annotation helpers
# ---------------------------------------------------------------------------

"Dashed guide line + label for a horizontal reference value (a forcing bound)."
function guide_h!(ax, y, label; color = GREY, xr, tx = nothing)
    hlines!(ax, y; color = color, linestyle = :dash, linewidth = 1.2)
    text!(ax, isnothing(tx) ? xr[1] : tx, y; text = label, color = color,
          align = (:left, :bottom), fontsize = 13, offset = (4, 2))
end

"Dashed vertical guide + label at an event time."
function guide_v!(ax, x, label, yr; color = GREY, align = (:left, :top))
    vlines!(ax, x; color = color, linestyle = :dash, linewidth = 1.2)
    text!(ax, x, yr[2]; text = label, color = color, align = align,
          fontsize = 13, offset = (4, -4))
end

"Double-tick horizontal span with a centered label (for dt_ramp, dt_conv, period)."
function span_h!(ax, x1, x2, y, label; color = ORANGE, tick = 0.25)
    lines!(ax, [x1, x2], [y, y]; color = color, linewidth = 2)
    lines!(ax, [x1, x1], [y - tick, y + tick]; color = color, linewidth = 2)
    lines!(ax, [x2, x2], [y - tick, y + tick]; color = color, linewidth = 2)
    text!(ax, 0.5 * (x1 + x2), y; text = label, color = color,
          align = (:center, :bottom), fontsize = 13, offset = (0, 4))
end

save_fig(name, fig) = (save(joinpath(OUT, name), fig); println("  wrote ", name))

# ---------------------------------------------------------------------------
# 1. Overview: all four time-driven methods
# ---------------------------------------------------------------------------

function fig_overview()
    t = range(0, 3000; length = 601)
    fig = Figure(size = (900, 560))

    ax1 = Axis(fig[1, 1]; title = "ramp-time", ylabel = "forcing  f")
    lines!(ax1, t, f_ramp_time.(t; dt_init = 0, dt_ramp = 1500, fmin = 0, fmax = 10);
           color = BLUE, linewidth = 2.5)

    ax2 = Axis(fig[1, 2]; title = "ramp-slope")
    lines!(ax2, t, f_ramp_slope.(t; dt_init = 0, df_dt_max = 0.0045, fmin = 0, fmax = 10);
           color = BLUE, linewidth = 2.5)

    ax3 = Axis(fig[2, 1]; title = "ramp-time-step", xlabel = "time  [yr]", ylabel = "forcing  f")
    lines!(ax3, t, f_ramp_time_step.(t; dt_init = 0, dt_ramp = 1000, dt_conv = 1000,
                                     fmin = 0, fmax = 10, fconv = 5); color = BLUE, linewidth = 2.5)

    ax4 = Axis(fig[2, 2]; title = "sin", xlabel = "time  [yr]")
    lines!(ax4, t, f_sin.(t; dt_ramp = 1500, fmin = 0, fmax = 10); color = BLUE, linewidth = 2.5)

    linkyaxes!(ax1, ax2, ax3, ax4)
    Label(fig[0, :], "tsgen — time-driven forcing methods";
          fontsize = 18, font = :bold, padding = (0, 0, 6, 0))
    save_fig("tsgen_overview.png", fig)
end

# ---------------------------------------------------------------------------
# 2. ramp-time (annotated)
# ---------------------------------------------------------------------------

function fig_ramp_time()
    dt_init, dt_ramp, fmin, fmax = 200.0, 1000.0, 0.0, 10.0
    t = range(0, 2000; length = 401)
    f = f_ramp_time.(t; dt_init, dt_ramp, fmin, fmax)
    xr, yr = (0.0, 2000.0), (-1.5, 12.0)

    fig = Figure(size = (820, 480))
    ax = Axis(fig[1, 1]; xlabel = "time  [yr]", ylabel = "forcing  f",
              title = "method = \"ramp-time\"", limits = (xr, yr))

    guide_h!(ax, fmax, "f_max"; xr, tx = 1450, color = RED)
    guide_h!(ax, fmin, "f_min"; xr, tx = 40,  color = RED)
    guide_v!(ax, dt_init, "dt_init\n(hold)", yr; color = GREEN)
    guide_v!(ax, dt_init + dt_ramp, "ramp ends", yr; color = GREEN)
    span_h!(ax, dt_init, dt_init + dt_ramp, 5.6, "dt_ramp")

    lines!(ax, t, f; color = BLUE, linewidth = 3)
    scatter!(ax, [dt_init, dt_init + dt_ramp], [fmin, fmax]; color = BLUE, markersize = 11)
    text!(ax, 720, 2.2; text = "rate = |f_max − f_min| / dt_ramp",
          color = GREY, fontsize = 13, align = (:left, :center))

    save_fig("tsgen_ramp_time.png", fig)
end

# ---------------------------------------------------------------------------
# 3. ramp-slope / const (annotated)
# ---------------------------------------------------------------------------

function fig_ramp_slope()
    dt_init, df_dt_max, fmin, fmax = 200.0, 0.008, 0.0, 10.0
    t_sat = dt_init + (fmax - fmin) / df_dt_max
    t = range(0, 2000; length = 401)
    f = f_ramp_slope.(t; dt_init, df_dt_max, fmin, fmax)
    xr, yr = (0.0, 2000.0), (-1.5, 12.0)

    fig = Figure(size = (820, 480))
    ax = Axis(fig[1, 1]; xlabel = "time  [yr]", ylabel = "forcing  f",
              title = "method = \"ramp-slope\"  (constant rate df_dt_max)",
              limits = (xr, yr))

    guide_h!(ax, fmax, "f_max  (clamp)"; xr, tx = 1350, color = RED)
    guide_v!(ax, dt_init, "dt_init", yr; color = GREEN)
    guide_v!(ax, t_sat, "saturates", yr; color = GREEN, align = (:right, :top))

    lines!(ax, t, f; color = BLUE, linewidth = 3)
    # slope triangle
    x0, x1 = 500.0, 900.0
    y0 = f_ramp_slope(x0; dt_init, df_dt_max, fmin, fmax)
    y1 = f_ramp_slope(x1; dt_init, df_dt_max, fmin, fmax)
    lines!(ax, [x0, x1, x1], [y0, y0, y1]; color = ORANGE, linewidth = 2)
    text!(ax, x1 + 30, 0.5 * (y0 + y1); text = "df_dt_max", color = ORANGE,
          fontsize = 13, align = (:left, :center))

    save_fig("tsgen_ramp_slope.png", fig)
end

# ---------------------------------------------------------------------------
# 4. ramp-time-step / triangle (annotated)
# ---------------------------------------------------------------------------

function fig_ramp_time_step()
    dt_init, dt_ramp, dt_conv = 0.0, 1000.0, 1000.0
    fmin, fmax, fconv = 0.0, 10.0, 5.0
    t = range(0, 3000; length = 601)
    f = f_ramp_time_step.(t; dt_init, dt_ramp, dt_conv, fmin, fmax, fconv)
    xr, yr = (0.0, 3000.0), (-1.5, 12.5)

    fig = Figure(size = (820, 480))
    ax = Axis(fig[1, 1]; xlabel = "time  [yr]", ylabel = "forcing  f",
              title = "method = \"ramp-time-step\"  (triangle)", limits = (xr, yr))

    guide_h!(ax, fmax,  "f_max"; xr, tx = 40, color = RED)
    guide_h!(ax, fconv, "f_conv"; xr, tx = 2450, color = RED)
    guide_v!(ax, dt_ramp, "peak", yr; color = GREEN)
    guide_v!(ax, dt_ramp + dt_conv, "hold", yr; color = GREEN)
    span_h!(ax, 0, dt_ramp, 11.4, "dt_ramp")
    span_h!(ax, dt_ramp, dt_ramp + dt_conv, 11.4, "dt_conv")

    lines!(ax, t, f; color = BLUE, linewidth = 3)
    scatter!(ax, [dt_ramp, dt_ramp + dt_conv], [fmax, fconv]; color = BLUE, markersize = 11)

    save_fig("tsgen_ramp_time_step.png", fig)
end

# ---------------------------------------------------------------------------
# 5. sin (annotated)
# ---------------------------------------------------------------------------

function fig_sin()
    dt_ramp, fmin, fmax = 2000.0, 0.0, 10.0   # dt_ramp is the period
    amp = 0.5 * (fmax - fmin)
    off = 0.5 * (fmax + fmin)
    t = range(0, 4000; length = 801)
    f = f_sin.(t; dt_ramp, fmin, fmax)
    xr, yr = (0.0, 4000.0), (-2.0, 12.5)

    fig = Figure(size = (820, 480))
    ax = Axis(fig[1, 1]; xlabel = "time  [yr]", ylabel = "forcing  f",
              title = "method = \"sin\"", limits = (xr, yr))

    guide_h!(ax, fmax, "f_max"; xr, tx = 40, color = RED)
    guide_h!(ax, fmin, "f_min"; xr, tx = 40, color = RED)
    guide_h!(ax, off,  "mean = ½(f_max+f_min)"; xr, tx = 2650, color = GREY)
    span_h!(ax, 500, 500 + dt_ramp, 11.6, "period = dt_ramp")
    # amplitude arrow
    lines!(ax, [500, 500], [off, fmax]; color = ORANGE, linewidth = 2)
    text!(ax, 520, 0.5 * (off + fmax); text = "amp = ½(f_max−f_min)",
          color = ORANGE, fontsize = 13, align = (:left, :center))

    lines!(ax, t, f; color = BLUE, linewidth = 3)
    save_fig("tsgen_sin.png", fig)
end

# ---------------------------------------------------------------------------
# 6. noise on top of the mean forcing (sigma)
# ---------------------------------------------------------------------------

function fig_noise()
    dt_init, dt_ramp, fmin, fmax = 0.0, 1500.0, 0.0, 10.0
    σ = 0.6
    t = range(0, 2000; length = 201)
    fmean = f_ramp_time.(t; dt_init, dt_ramp, fmin, fmax)
    # deterministic pseudo-noise (reproducible, no RNG dependency)
    noise = σ .* sin.(0.06 .* t) .* cos.(0.017 .* t .+ 1.3)
    fnow = fmean .+ noise

    fig = Figure(size = (820, 460))
    ax = Axis(fig[1, 1]; xlabel = "time  [yr]", ylabel = "forcing  f",
              title = "noise: f_now = f_mean + eps,   eps ~ N(0, sigma)",
              limits = ((0.0, 2000.0), (-1.5, 12.0)))

    band!(ax, t, fmean .- σ, fmean .+ σ; color = (BLUE, 0.12))
    lines!(ax, t, fnow;  color = (GREY, 0.9), linewidth = 1.2, label = "f_now  (mean + noise)")
    lines!(ax, t, fmean; color = BLUE, linewidth = 3, label = "f_mean")
    text!(ax, 1550, 8.4; text = "±sigma", color = BLUE, fontsize = 13, align = (:left, :center))
    axislegend(ax; position = :lt, framevisible = false)

    save_fig("tsgen_noise.png", fig)
end

# ---------------------------------------------------------------------------
# 7. series: clamped-linear interpolation of a tabulated curve
# ---------------------------------------------------------------------------

function fig_series()
    tt = [0.0, 1000.0, 2000.0]          # table times  (matches test_tsgen series)
    vv = [0.0, 10.0, 30.0]              # table values
    xr = (-400.0, 2600.0)

    # clamped-linear interpolation
    function interp(x)
        x <= tt[1]   && return vv[1]
        x >= tt[end] && return vv[end]
        i = searchsortedlast(tt, x)
        w = (x - tt[i]) / (tt[i+1] - tt[i])
        vv[i] + w * (vv[i+1] - vv[i])
    end
    xs = range(xr...; length = 601)

    fig = Figure(size = (820, 460))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "value",
              title = "series — clamped-linear interpolation of a tabulated curve",
              limits = (xr, (-4.0, 34.0)))

    # clamped regions (held flat outside the table range)
    vlines!(ax, [tt[1], tt[end]]; color = GREY, linestyle = :dash, linewidth = 1.1)
    lines!(ax, [xr[1], tt[1]], [vv[1], vv[1]]; color = ORANGE, linewidth = 2.5)
    lines!(ax, [tt[end], xr[2]], [vv[end], vv[end]]; color = ORANGE, linewidth = 2.5)
    text!(ax, xr[1] + 30, vv[1] + 1.5; text = "clamped\n(held at endpoint)",
          color = ORANGE, fontsize = 12, align = (:left, :bottom))
    text!(ax, 0.5 * (tt[end] + xr[2]), vv[end]; text = "clamped",
          color = ORANGE, fontsize = 12, align = (:center, :bottom), offset = (0, 6))

    lines!(ax, xs, interp.(xs); color = BLUE, linewidth = 2.5)
    scatter!(ax, tt, vv; color = BLUE, markersize = 14)
    text!(ax, tt, vv; text = ["(0, 0)", "(1000, 10)", "(2000, 30)"],
          color = BLUE, fontsize = 12, align = (:left, :top), offset = (8, -8))

    save_fig("series_interp.png", fig)
end

# ---------------------------------------------------------------------------
# 8. root_finder: Newton's method
# ---------------------------------------------------------------------------

function fig_newton()
    f(x)  = x^2 - 2.0            # root at sqrt(2)
    fp(x) = 2.0x
    x = 2.6                       # initial guess
    iters = [x]
    for _ in 1:3
        x = x - f(x) / fp(x)
        push!(iters, x)
    end

    xs = range(0.6, 2.9; length = 400)
    fig = Figure(size = (760, 480))
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "f(x)",
              title = "Newton's method:  x ← x − f(x) / f′(x)",
              limits = ((0.6, 2.9), (-3.0, 6.5)))

    hlines!(ax, 0; color = GREY, linewidth = 1)
    lines!(ax, xs, f.(xs); color = BLUE, linewidth = 2.5, label = "f(x) = x² − 2")

    for i in 1:3
        x0 = iters[i]; y0 = f(x0); x1 = iters[i+1]
        lines!(ax, [x0, x0], [0, y0]; color = (GREY, 0.6), linestyle = :dot, linewidth = 1)
        lines!(ax, [x0, x1], [y0, 0]; color = ORANGE, linewidth = 1.6)   # tangent to x-intercept
        scatter!(ax, [x0], [y0]; color = BLUE, markersize = 9)
        text!(ax, x0, -0.5; text = "x$(i-1)", color = GREY, fontsize = 12, align = (:center, :top))
    end
    scatter!(ax, [sqrt(2)], [0]; color = RED, markersize = 12)
    text!(ax, sqrt(2), 0.35; text = "root = √2", color = RED, fontsize = 12, align = (:center, :bottom))
    axislegend(ax; position = :lt, framevisible = false)

    save_fig("root_finder_newton.png", fig)
end

# ---------------------------------------------------------------------------
# 9. staggering: Arakawa C-grid
# ---------------------------------------------------------------------------

function fig_cgrid()
    fig = Figure(size = (620, 560))
    ax = Axis(fig[1, 1]; title = "Arakawa C-grid staggering",
              aspect = DataAspect(), limits = ((-0.6, 3.6), (-0.6, 3.6)))
    hidedecorations!(ax); hidespines!(ax)

    n = 3
    # cell grid lines
    for k in 0:n
        lines!(ax, [0, n], [k, k]; color = (GREY, 0.4), linewidth = 1)
        lines!(ax, [k, k], [0, n]; color = (GREY, 0.4), linewidth = 1)
    end

    centers = [(i + 0.5, j + 0.5) for i in 0:n-1, j in 0:n-1]
    uxs = [(i, j + 0.5) for i in 0:n, j in 0:n-1]        # x-edges  (u / ac-x)
    vys = [(i + 0.5, j) for i in 0:n-1, j in 0:n]        # y-edges  (v / ac-y)

    sc_aa = scatter!(ax, first.(centers)[:], last.(centers)[:];
                     color = BLUE, markersize = 15)
    sc_u = scatter!(ax, first.(uxs)[:], last.(uxs)[:];
                    color = ORANGE, marker = :rtriangle, markersize = 15)
    sc_v = scatter!(ax, first.(vys)[:], last.(vys)[:];
                    color = GREEN, marker = :utriangle, markersize = 15)

    Legend(fig[2, 1],
        [sc_aa, sc_u, sc_v],
        ["aa-node — scalar (H, T, …)", "ac-node x — u component", "ac-node y — v component"],
        orientation = :horizontal, framevisible = false, nbanks = 3)

    save_fig("staggering_cgrid.png", fig)
end

# ---------------------------------------------------------------------------

function main()
    println("Generating fesm-utils documentation figures ...")
    fig_overview()
    fig_ramp_time()
    fig_ramp_slope()
    fig_ramp_time_step()
    fig_sin()
    fig_noise()
    fig_series()
    fig_newton()
    fig_cgrid()
    println("Done.")
end

main()
