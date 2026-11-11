# Light-Written Photonics Lab

Open `index.html` in a browser (single file, no build). First load takes about 15 s while it simulates 26 MZIs.

**Idea:** instead of etching waveguides, a Fresnel plate (phase hologram) projects pump light onto a photosensitive slab. Where the light lands, the refractive index changes, and that index pattern becomes the optical circuit. Rewriting the plate rewrites the processor.

## What it computes

1. **Plate pattern:** Gerchberg–Saxton (noise-region variant) computes the phase pattern for the plate, including a lens term so no lens is needed. It then measures the projected light's speckle, background and efficiency, with M frames averaged.
2. **Written index:** n(x,z) = n₀ + Δn·pump, blurred by the projector PSF (0.22 λₚ/NA) and material diffusion, with the measured speckle and a slow illumination drift added on top.
3. **Components, by split-step beam propagation (BPM) at 1550 nm:**
   - Waveguides designed single-mode (V = 0.8·π/2).
   - Directional couplers set to 50:50 using supermodes, then trimmed with BPM including the S-bends.
   - Mach–Zehnder switches: extinction ratio, insertion loss, and the drive needed for π (electro-optic V, thermo-optic ΔT, or extra pump light).
   - Ring add-drop filters, with coupling from coupled-mode theory along the curved gap and bend loss from BPM in a conformally mapped bend.
4. **Processor:** a 2-bit ripple-carry adder in directed logic (11 MZIs, dual-rail carry). Every switch is its own BPM run with its own noise. All 16 inputs are checked, and the eye opening is reported for each output.

Each step runs for a light-written slab (SBN, LiNbO₃:Fe, photopolymer, or optically pumped silicon) and for a fabricated silica PLC reference (Δn = 7.4×10⁻³).

## Results with default settings (SBN, Δn = 10⁻³)

| | Light-written | Fabricated silica |
|---|---|---|
| MZI length | 26.9 mm | 5.6 mm |
| Extinction (bar / cross) | 40 / 36 dB | 43 / 40 dB |
| Insertion loss per MZI | 1.4 dB (material loss) | 0.2 dB |
| Drive for π | 5.5 V (electro-optic) | ΔT ≈ 43 K |
| 2-bit adder | 16/16, eye ≥ 27 dB | 16/16, eye ≥ 33 dB |
| Ring filter | unusable: R ≈ 170 mm needed | R ≈ 5.5 mm, Q ≈ 10⁵ |

Weak light-induced Δn works for interferometric logic (MZIs). Rings need Δn ≳ 10⁻², which points to phase-change films.

Physics notes, limitations and sources are in the page's tab 4.
