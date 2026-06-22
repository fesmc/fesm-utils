# Julia example to perform spherical harmonic transforms using SHTns library
# install SHTns with:
#    import Pkg; Pkg.add("SHTns")
import SHTns

lmax = 13   # maximum degree
sh = SHTns.SHTnsCfg(13,13,1)
theta, phi = SHTns.grid(sh)  # get physical grid

alm = zeros(Complex{Float64}, sh.nlm)   # an array that can hold the spherical harmonic coefficients
lm10 = SHTns.LM(sh, 1,0)  # index in array for l=1, m=0
alm[lm10] = sqrt(4*pi/3)  # l=1,m=0

# synthesis: from spectral qlm to spatial x
x=SHTns.synth(sh, alm)

# as we used l=1,m=0, it should be the same as cos(theta) (sh.ct):
display( maximum(abs.(x[:,1] - sh.ct)) < 1e-15 )

# analysis: transform back spatial data to spectral
alm2 = SHTns.analys(sh, x)
test = isapprox(alm2,alm)  # alm and alm2 should match, within machine precision
if test
    println("OK")
end

import Plots
# plot an iso-longitude cut (meridional) of phi=0
Plots.plot(acos.(sh.ct), x[:,1])  # should be cos(theta) between 0 and pi

