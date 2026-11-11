# Steering light with a second wave

An interactive simulation, and the physics behind it, for one question: Feynman explains refraction as the incident wave plus waves re-radiated by the material. Can you get the same bending by projecting a wave of your own onto the light?

Open `index.html` in a browser. It is a single file with no build step.

## Short answer

| Way of adding the wave | Does the beam change direction? | Why |
|---|---|---|
| A second free beam crossing it in vacuum | **No** | Maxwell's equations are linear in vacuum. The beams interfere where they overlap and then pass through each other unchanged. |
| A **phase-locked** sheet of sources, same frequency, with a phase gradient across it | **Yes** | This is what the material does (Feynman Ch. 31). It obeys sin θₜ = sin θᵢ + (λ/2π)·dΦ/dx. |
| The same sheet with a frequency offset Δf | **No** (on average) | The interference term averages to zero over a time 1/Δf. |
| A wave that **writes a grating** into a medium (sound wave, or two pump beams) | **Yes** | The grating equation sin θₘ = sin θᵢ + mλ/Λ. Acousto-optic deflectors and photorefractive / Kerr gratings work this way. |

## The governing relationship

Along the sheet, the transmitted wave has to match the phase of the total field:

```
k sin θt = k sin θi + dΦ/dx + m·2π/Λ
sin θt   = sin θi + (λ/2π)·dΦ/dx + m·λ/Λ
```

- Linear ramp Φ = (Δφ/d)·x: sin θₜ = sin θᵢ + λΔφ/(2πd). This is the generalized Snell's law and the phased-array scan equation. Grating lobes appear at + mλ/d when d > λ/(1+|sin θ|).
- Periodic Φ = M sin(2πx/Λ): sin θₘ = sin θᵢ + mλ/Λ, with order powers J_m(M)² for a thin grating. Λ = v/f for a sound wave, and Λ = λₚ/(2 sin α) for two pump beams.
- Uniform Φ: no turning, only a phase delay δ. This equals a slab with (n−1)Δz = δλ/2π, which is Feynman's plate.

The full write-up, the model details and the sources are in the "notes" section of `index.html`.

## Sources

1. Feynman, Leighton, Sands, *The Feynman Lectures on Physics* Vol. I, Ch. 31 (refractive index) and Ch. 30 (§30-2 grating, §30-7 field of a plane of charges). https://www.feynmanlectures.caltech.edu/I_31.html
2. Born & Wolf, *Principles of Optics*, Ch. 2 (§2.4): Ewald–Oseen extinction theorem.
3. N. Yu et al., "Light propagation with phase discontinuities: generalized laws of reflection and refraction," *Science* 334, 333 (2011). doi:10.1126/science.1210713
4. C. Pfeiffer & A. Grbic, "Metamaterial Huygens' surfaces," *PRL* 110, 197401 (2013). doi:10.1103/PhysRevLett.110.197401
5. R. J. Mailloux, *Phased Array Antenna Handbook* (scan equation, grating lobes).
6. Saleh & Teich, *Fundamentals of Photonics*, acousto-optics chapter (Bragg condition, Λ = v/f).
7. C. V. Raman & N. S. N. Nath, *Proc. Indian Acad. Sci. A* 2, 406 (1935): Bessel order intensities.
8. H. Kogelnik, *Bell Syst. Tech. J.* 48, 2909 (1969): thick-grating Bragg selectivity.
9. P. Yeh, *Introduction to Photorefractive Nonlinear Optics* (1993): light-written gratings.
10. Heisenberg & Euler, *Z. Phys.* 98, 714 (1936); ATLAS, *Nature Physics* 13, 852 (2017): vacuum light-by-light scattering, negligible at laser intensities.
