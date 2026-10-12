# zremap_figures.jl
# ---------------------------------------------------------------------------
# Figures of the zremap documentation page, from the output of the demo
# program (present-day Greenland SMB remapped to the surfaces of the ICE-6G_C
# deglaciation):
#
#   cd fesm-utils
#   make demo-zremap
#   ./bin/demo_zremap_grl.x docs/figures/zremap/zremap_grl.nml
#   cd docs/figures
#   julia --project=. zremap_figures.jl
#
# Outputs are committed PNGs so the site builds without Julia or the data.
# ---------------------------------------------------------------------------

using CairoMakie
using NCDatasets
CairoMakie.activate!(px_per_unit = 2)

const OUT  = @__DIR__
const FILE = joinpath(OUT, "zremap", "zremap_grl.nc")

const BLUE   = colorant"#1f6feb"
const ORANGE = colorant"#e8710a"
const GREY   = colorant"#57606a"
const RED    = colorant"#d1005c"

set_theme!(Theme(fontsize = 14))

save_fig(name, fig) = (save(joinpath(OUT, name), fig); println("  wrote ", name))

# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

function load(file)
    ds = NCDataset(file)
    get(v)  = Float64.(coalesce.(ds[v][:, :], NaN))
    get3(v) = Float64.(coalesce.(ds[v][:, :, :], NaN))
    d = (
        xc = Float64.(ds["xc"][:]), yc = Float64.(ds["yc"][:]),
        zb = Float64.(ds["zb"][:]), time = Float64.(ds["time"][:]),
        f_ref = get("f_ref"), z_ref = get("z_ref"),
        w = get("w"), dfdz = get("dfdz"),
        mask_ref = ds["mask_ref"][:, :] .== 1,
        f_tab = get3("f_tab"), z_lo = get("z_lo"), z_hi = get("z_hi"),
        # per time slice (x, y, time)
        z_now = get3("z_now"), f_now = get3("f_now"), f_grad = get3("f_grad"),
        dz_ext = get3("dz_ext"), mask_now = ds["mask_now"][:, :, :] .== 1,
    )
    close(ds)
    d
end

"The fields of time slice t [kyr] (2D, for the map and profile figures)."
function slice(d, t)
    k = argmin(abs.(d.time .- t))
    merge(d, (z_now = d.z_now[:, :, k], f_now = d.f_now[:, :, k], f_grad = d.f_grad[:, :, k],
              dz_ext = d.dz_ext[:, :, k], mask_now = d.mask_now[:, :, k]))
end

const SMBCOL = (colormap = :RdBu, colorrange = (-2000, 2000))

masked(v, m) = ifelse.(m, v, NaN)

function map_axis(f, title)
    ax = Axis(f; title = title, aspect = DataAspect(), titlesize = 14)
    hidedecorations!(ax)
    hidespines!(ax)
    ax
end

outline!(ax, d, m; color = :black, linewidth = 0.8, linestyle = :solid) =
    contour!(ax, d.xc, d.yc, Float64.(m); levels = [0.5], color = color,
             linewidth = linewidth, linestyle = linestyle)

# ---------------------------------------------------------------------------
# 1. Maps: reference, LGM remap, constant-gradient remap; support and dz_ext
# ---------------------------------------------------------------------------

function fig_maps(d; pts = [])
    fig = Figure(size = (1000, 900))
    smbcol = SMBCOL
    land = .!isnan.(d.f_ref)

    ax = map_axis(fig[1, 1], "(a) Present-day SMB (MAR)")
    heatmap!(ax, d.xc, d.yc, masked(d.f_ref, d.mask_ref); smbcol...)
    outline!(ax, d, land; color = GREY, linewidth = 0.5)
    for (k, (i, j)) in enumerate(pts)
        scatter!(ax, [d.xc[i]], [d.yc[j]]; color = :black, markersize = 7)
        text!(ax, d.xc[i], d.yc[j]; text = string('A' + k - 1), offset = (5, 0),
              align = (:left, :center), fontsize = 14, font = :bold)
    end

    ax = map_axis(fig[1, 2], "(b) zremap SMB on the LGM surface")
    hm = heatmap!(ax, d.xc, d.yc, masked(d.f_now, d.mask_now); smbcol...)
    outline!(ax, d, d.mask_ref)
    outline!(ax, d, land; color = GREY, linewidth = 0.5)

    ax = map_axis(fig[1, 3], "(c) Constant gradient dfdz_lo")
    heatmap!(ax, d.xc, d.yc, masked(d.f_grad, d.mask_now); smbcol...)
    outline!(ax, d, d.mask_ref)
    outline!(ax, d, land; color = GREY, linewidth = 0.5)
    Colorbar(fig[1, 4], hm; label = "SMB  [mm w.e. a⁻¹]")

    ax = map_axis(fig[2, 1], "(d) LGM − present surface")
    hm = heatmap!(ax, d.xc, d.yc, masked(d.z_now .- d.z_ref, d.mask_now);
                  colormap = :viridis, colorrange = (0, 1500))
    outline!(ax, d, d.mask_ref; color = :white)
    outline!(ax, d, land; color = GREY, linewidth = 0.5)
    Colorbar(fig[3, 1], hm; vertical = false, flipaxis = false, label = "Δz  [m]",
             ticks = 0:500:1000)

    ax = map_axis(fig[2, 2], "(e) Support W")
    hm = heatmap!(ax, d.xc, d.yc, d.w; colormap = :magma, colorrange = (0, 1))
    contour!(ax, d.xc, d.yc, d.w; levels = [0.25], color = :white,
             linestyle = :dash, linewidth = 1.2)
    outline!(ax, d, d.mask_now; color = BLUE, linewidth = 1.2)
    outline!(ax, d, land; color = GREY, linewidth = 0.5)
    Colorbar(fig[3, 2], hm; vertical = false, flipaxis = false,
             label = "W  [1]   (dashed: W = w_full)", ticks = 0.25:0.25:0.75)

    ax = map_axis(fig[2, 3], "(f) Elevation beyond the data, dz_ext")
    hm = heatmap!(ax, d.xc, d.yc, masked(d.dz_ext, d.mask_now);
                  colormap = Reverse(:PuOr), colorrange = (-1000, 1000))
    outline!(ax, d, d.mask_ref)
    outline!(ax, d, land; color = GREY, linewidth = 0.5)
    Colorbar(fig[3, 3], hm; vertical = false, flipaxis = false, label = "dz_ext  [m]",
             ticks = -500:500:500)

    rowgap!(fig.layout, 2, 4)
    fig
end

# ---------------------------------------------------------------------------
# 2. Profiles at a few points: data, profile, LGM elevation
# ---------------------------------------------------------------------------

function fig_profiles(d, pts; sigma = 100.0)
    fig = Figure(size = (1000, 380))
    for (k, (i, j)) in enumerate(pts)
        ax = Axis(fig[1, k]; title = "$(('A' + k - 1))   W = $(round(d.w[i, j], digits = 2))",
                  xlabel = "SMB  [mm w.e. a⁻¹]", ylabel = k == 1 ? "elevation  [m]" : "")

        # data around the point, shaded by kernel weight
        r2 = ((d.xc .- d.xc[i]) .^ 2 .+ (d.yc' .- d.yc[j]) .^ 2) ./ sigma^2
        wk = exp.(-0.5 .* r2)
        sel = d.mask_ref .& (r2 .<= 9)
        scatter!(ax, d.f_ref[sel], d.z_ref[sel]; color = wk[sel], colormap = [:white, GREY],
                 colorrange = (0, 1), markersize = 5)

        # profile: within the data range solid, extrapolated dashed
        p   = d.f_tab[i, j, :]
        lo, hi = d.z_lo[i, j], d.z_hi[i, j]
        if !isnan(lo)
            inr = (d.zb .>= lo) .& (d.zb .<= hi)
            lines!(ax, p, d.zb; color = BLUE, linewidth = 1.5, linestyle = :dash)
            lines!(ax, p[inr], d.zb[inr]; color = BLUE, linewidth = 2.5)
        end

        # reference and LGM elevation at the point
        hlines!(ax, d.z_ref[i, j]; color = GREY, linestyle = :dot)
        hlines!(ax, d.z_now[i, j]; color = ORANGE, linewidth = 1.5)
        scatter!(ax, [d.f_now[i, j]], [d.z_now[i, j]]; color = ORANGE, markersize = 11)
        scatter!(ax, [d.f_grad[i, j]], [d.z_now[i, j]]; color = RED, marker = :xcross,
                 markersize = 11)

        ax.xticks = -4000:2000:2000
        xlims!(ax, -4000, 2000)
        ylims!(ax, 0, 4000)
    end
    Legend(fig[1, length(pts) + 1],
           [MarkerElement(color = GREY, marker = :circle),
            LineElement(color = BLUE, linewidth = 2.5),
            LineElement(color = BLUE, linestyle = :dash),
            LineElement(color = GREY, linestyle = :dot),
            MarkerElement(color = ORANGE, marker = :circle),
            MarkerElement(color = RED, marker = :xcross)],
           ["data (shade: kernel weight)", "profile (data range)", "profile (extrapolated)",
            "z_ref", "zremap at z_now (LGM)", "constant gradient"];
           framevisible = false, labelsize = 12)
    fig
end

# ---------------------------------------------------------------------------
# 3. Deglaciation: SMB on the ice of a few time slices, and integrated SMB
# ---------------------------------------------------------------------------

function fig_progression(d; times = [-21, -16, -13, -11, -9, 0], dx = 16e3)
    fig = Figure(size = (1000, 540))
    land = .!isnan.(d.f_ref)
    hm = nothing
    for (n, t) in enumerate(times)
        s = slice(d, t)
        ax = map_axis(fig[1, n], "$(Int(-t)) ka")
        hm = heatmap!(ax, d.xc, d.yc, masked(s.f_now, s.mask_now); SMBCOL...)
        outline!(ax, d, d.mask_ref; linewidth = 0.6)
        outline!(ax, d, land; color = GREY, linewidth = 0.5)
    end
    Colorbar(fig[1, length(times) + 1], hm; label = "SMB  [mm w.e. a⁻¹]", height = Relative(0.85))

    # Integrated SMB over the ice of each slice [Gt a-1] (cells of dx^2)
    gt(f, m) = sum(filter(!isnan, f[m])) * 1e-3 * dx^2 * 1e3 / 1e12
    nt = length(d.time)
    smb_z = [gt(d.f_now[:, :, k], d.mask_now[:, :, k]) for k in 1:nt]
    smb_g = [gt(d.f_grad[:, :, k], d.mask_now[:, :, k]) for k in 1:nt]
    area  = [count(d.mask_now[:, :, k]) * dx^2 / 1e12 for k in 1:nt]

    ga = fig[2, 1:length(times) + 1] = GridLayout()
    ax = Axis(ga[1, 1]; xlabel = "time  [ka]", ylabel = "SMB  [Gt a⁻¹]",
              xreversed = true)
    lines!(ax, -d.time, smb_z; color = BLUE, linewidth = 2.5, label = "zremap")
    lines!(ax, -d.time, smb_g; color = RED, linewidth = 2, linestyle = :dash,
           label = "constant gradient dfdz_lo")
    hlines!(ax, gt(d.f_ref, d.mask_ref); color = GREY, linestyle = :dot,
            label = "present-day MAR")
    vlines!(ax, -Float64.(times); color = (GREY, 0.3), linewidth = 1)
    axislegend(ax; position = :lt, framevisible = false, labelsize = 12)

    ax2 = Axis(ga[1, 2]; xlabel = "time  [ka]", ylabel = "ice area  [10⁶ km²]",
               xreversed = true)
    lines!(ax2, -d.time, area; color = GREY, linewidth = 2)
    colsize!(ga, 2, Relative(0.3))
    rowsize!(fig.layout, 2, Fixed(210))
    fig
end

# ---------------------------------------------------------------------------

d   = load(FILE)
lgm = slice(d, -21)

# Profile points (grid indices): A south-west ablation zone, B north-west
# margin, C the LGM ice point beyond the present ice that lies furthest above
# the data range with some support (W > 0.1)
cand = lgm.mask_now .& .!d.mask_ref .& (d.w .> 0.1) .& .!isnan.(lgm.dz_ext)
kc   = argmax(ifelse.(cand, lgm.dz_ext, -Inf))
const PTS = [(30, 50), (33, 150), (kc[1], kc[2])]

save_fig("zremap_maps.png", fig_maps(lgm; pts = PTS))
save_fig("zremap_profiles.png", fig_profiles(lgm, PTS))
save_fig("zremap_progression.png", fig_progression(d))
