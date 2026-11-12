# Optical GPU, written with light

A design tool and compiler for **optical tensor processors whose circuits and weights are written by projected light** into a photosensitive film. This is v3 of the light-steering study: `light-steering-sim/` (v1) → `photonic-lab/` (v2) → here.

Each tensor core computes y = W·x for an arbitrary N×N real matrix in one pass of light: a Vᵀ mesh, then Σ attenuators, then a U mesh of MZIs (SVD + Clements decomposition). An "optical GPU" is an array of such cores driven by an electronic host. The default design `OGPU-1` has four 16×16 cores on a Sb₂Se₃-on-SiN film.

## Use it

| What | How |
|---|---|
| Design in the browser | Open `designer/index.html` (single file). Edit blocks, cores, platform and process, then press *Compile & simulate*. |
| Compile a design file into manufacturing outputs | `node tools/run.js designs/ogpu-1.json build/ogpu-1 --all-tiles` |
| Rebuild the designer after editing the engine | `node tools/bundle.js` |
| Build the hardware | `docs/BUILD.md` |

No dependencies beyond Node ≥ 18 and a browser.

## What the compiler produces (`build/ogpu-1/`, committed as run)

| File | Content |
|---|---|
| `report.md` | Cell geometry, accuracy per core, performance vs silicon photonics and H100, die and writing time |
| `design.resolved.json` | Design + derived geometry, measured device pool, accuracy, performance |
| `layout.gds` | GDSII. Layer 1 waveguides, 2 weight patches (datatype = dose), 3 facet converters / terminators, 10 fiducials, 63 outline. Opens in KLayout |
| `program/core*_weights.csv` | Weight-patch list for the write–verify loop |
| `plate/dmd/` | **Programmable plate:** one 8-bit DMD frame per exposure tile + stage map `tiles.json` |
| `plate/lcos/` | **Fresnel plates:** phase holograms with a lens term, for a phase SLM or as an etched DOE; targets and simulated projections |
| `overview.png`, `zoom_mzi.png` | Raster previews |

## Engine (`core/otc.js`)

* **Device physics:** finite-difference mode solver; split-step BPM; conformal-map bend loss; 50:50 coupler calibration; process centering on written test MZIs.
* **Noisy writes:** each write includes blur (projector NA + material diffusion), speckle (plate type, frames averaged) and illumination drift. The coupler angles, arm phase and loss of every written MZI are extracted from these simulations.
* **Programming:** SVD (one-sided Jacobi) with degenerate-basis randomization. Clements decomposition, exact to 10⁻¹⁵. Hardware-aware per-MZI phase fitting with phase-debt propagation (calibration).
* **System:** quantized weights, coherent detection, receiver shot noise from the power budget, DAC/ADC quantization → ENOB. Power and TOPS/W model.
* **Compilation:** layout (raised-cosine bends at the bend-loss radius, fan-in/out to 127 µm fiber pitch, Σ spread with dump terminators, graded facet converters). Rasterizer → exposure tiles. Gerchberg–Saxton holograms → Fresnel plate. GDSII and PNG writers.

## Honest limits

2-D effective-index model; Sb₂Se₃ contrast and stability are literature estimates; precision is about 5 bits; weights update in ms–s. The host does memory, nonlinearity and control. See the Notes tab in the designer and `docs/BUILD.md` §9.
