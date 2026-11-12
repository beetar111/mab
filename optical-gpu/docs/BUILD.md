# Building an optically written tensor processor

This guide turns a compiled design (`build/<design>/`) into a working chip and system. Numbers in brackets are from the reference build `build/ogpu-1/report.md`; your design's report gives its own values.

The process is maskless. The only lithography-grade step is the light pattern, and the circuit, the weights and every correction to them are written by the same station.

> **Status.** Each step uses equipment and materials that exist today. Each one has also been shown separately in the literature: projected-light waveguides, laser-written rewritable phase-change waveguides, programmable MZI meshes with hardware error correction, and coherent optical matrix multiplication. Combining all of them in one die is the open engineering task. Treat the process windows below as starting points to calibrate, not as qualified recipes.

---

## 0. What you are building

```
1550 nm laser ──► splitter ──► N×K modulators (DAC) ──► fiber array ──► ┌───────── written die ─────────┐ ──► fiber array ──► N×K coherent receivers (ADC) ──► host
                     │                                                    │ fan-in → Vᵀ mesh → Σ → U mesh → fan-out │                      ▲
                     └──────────────────────── local oscillator ──────────────────────────────────────────────────────┘
```

* **Die:** an unpatterned film stack. Waveguides, couplers, MZI meshes and weights exist only as written index changes [69.6 × 9.7 mm; 4 cores × 16×16; 1024 MZIs; 2048 weight patches].
* **Off-chip:** laser, modulators, coherent receivers, data converters and an FPGA host, all commercial parts. The die needs no metal and no electronics.

## 1. Blank wafers (standard foundry steps, no patterning)

Default platform `pcm_sin`:

| Layer | Process | Notes |
|---|---|---|
| Si wafer, 100 mm | n/a | Prime grade, < 1 µm bow |
| 4 µm SiO₂ | Wet thermal oxide | Lower cladding; keeps the mode off the Si |
| 300 nm Si₃N₄ | LPCVD, anneal 1100 °C | Unpatterned slab that carries the light vertically |
| 20 nm Sb₂Se₃ | RF sputter, room temperature (amorphous as deposited) | The writable layer. Crystalline areas have higher index |
| 10 nm SiO₂ | ALD or sputter | Cap: prevents oxidation and Se loss during writing |

Then:

1. Measure the slab's effective index and the film thickness on witness pieces (prism coupler or ellipsometry). The design assumes n_eff ≈ 1.80 and a full crystallization contrast Δn_eff ≈ 0.03. **Re-run the designer with your measured values.**
2. Dice the wafer into blanks 2 mm larger than the die, then polish both end facets (input/output) to optical quality. Facet angle tolerance is ±0.5°.

The write-once alternative (`polymer`) uses a glass wafer, a 5 µm low-index cladding and a 15 µm holographic photopolymer film, with no vacuum steps.

## 2. Writing station

| Item | Specification | Why |
|---|---|---|
| Writing laser | 532 nm, ≥ 1 W CW (crystallize), plus ns-pulsed 532 nm for erase | Sb₂Se₃ absorbs strongly in the visible |
| Programmable plate | **DMD** 2560×1600 class (8-bit gray by PWM) for production, **or** a phase-only **LCoS** SLM 4K class (3.74 µm pitch) for holographic writing | See "Which plate" below |
| Projection | 4f relay plus objective, NA ≥ 0.4, demagnification = plate pitch ÷ slab pixel [DMD ×11.3, slab pixel 0.667 µm = core width / 4] | The waveguide must be resolved with ≥ 4 pixels |
| Stage | XY air-bearing, ≥ 80 mm travel, ≤ 50 nm repeatability, laser-interferometer or encoder feedback; Z autofocus | Tile stitching error must be ≪ core width (2.7 µm) |
| Metrology in the station | Fiducial camera; power meter after the objective; a 1550 nm probe laser and an IR camera viewing the output facet | Closed-loop writing (step 5) |
| Environment | Class 1000 or better; temperature ±0.1 °C | Drift becomes coupler imbalance |

**Which plate.**

* **DMD (direct imaging)** is the production path. Illuminate it with an incoherent source (LED or speckle-reduced laser) to avoid speckle, which gives about 1% index noise. The frames are the PNGs in `plate/dmd/`, one per tile, with stage positions in `tiles.json`.
* **LCoS phase SLM** displays the Fresnel plates in `plate/lcos/`. It puts nearly all the laser power into the pattern and can focus into several planes, for future 3-D circuits. Average ≥ 8 frames, each computed from a different random start, to bring speckle from about 6% to about 2%.
* **Printed static plate:** etch the same phase pattern into fused silica as a 4-, 8- or 16-level DOE (grayscale or multi-mask lithography, from the hologram PNG or the GDS weight-free layer). One flash writes a whole tile, like the phase masks used to write fiber Bragg gratings. Use it for waveguides in volume, then trim the weights with the programmable plate.

## 3. Calibrate the process (once per film batch)

1. **Dose → index.** Write a test field of 16 stripes at doses from 0 to 100% and measure Δn with an interferometric microscope or written test MZIs. Fit a monotonic curve. Load its inverse as the gray-level lookup of the plate, so the design file's "dose level" means a known Δn.
2. **Line width and blur.** Write single lines and measure their profile (near-field imaging or AFM phase). Confirm the blur σ is ≤ the design value [0.33 µm].
3. **Coupler centering.** Write 10 test MZIs from the design, measure their cross-state extinction and fit the coupler angles. If the mean is not π/4, adjust the coupler length in the design (`OTC.designCell` does exactly this in simulation [+1 µm]) and recompile.
4. **Stitching.** Write a vernier pattern across tile boundaries and correct the stage or magnification until the overlay error is < 100 nm.

## 4. Write the die

1. Load the blank and level it with autofocus. Write the four fiducials from the GDS (layer 10) and image them to set the die coordinate frame.
2. Expose the tiles in serpentine order. Use a 2 µm overlap blended linearly across the seam to avoid a step in index. Waveguides are written at 50% of full contrast and weight patches at 50–100% (GDS layer 2, datatype = dose 0–255) [224 of the 410 DMD tiles contain waveguides, ≈3 min per die at 0.8 s per tile].
3. Write the facet mode converters last (GDS layer 3). The dose ramps from 12% at the facet to 50% over 100 µm, so the mode widens to match the fiber.

## 5. Package

1. Actively align two V-groove fiber arrays (127 µm pitch, PM fiber, [64 channels each]) to the input and output facets while monitoring transmission through a straight test lane. Cure with low-shrink UV epoxy.
2. Mount on a TEC with a thermistor, holding 25.00 ± 0.05 °C. Seal the lid with a nitrogen fill.
3. Expected insertion loss from fiber to fiber: 2 × facet loss (≈2 dB each) + material and device loss [≈9 dB through 33 MZI columns] [total ≈13 dB].

## 6. Calibrate the cores (write–verify): mandatory

The simulation shows that an uncalibrated mesh does not compute at all (ENOB < 0). Calibration is part of manufacturing:

1. **Characterize each MZI** in mesh order using the standard self-configuring sequence: earlier MZIs set to known states, inject on one port, read the output powers. Fit each MZI's two coupler angles and arm-phase offset [expected α 0.76–0.81 rad].
2. **Re-solve the phases for the real hardware** (`fitMZI`): each MZI gets (θ, φ) so that it matches the ideal block up to output phases. Those leftover phases are carried forward to later shifters and to the receiver phases.
3. **Write the weights closed-loop.** Expose each weight patch in small dose steps (crystallization only increases the index), measuring the MZI after each step, until it is within 2π/256 of target. If you overshoot, re-amorphize the patch with one ns pulse frame and repeat.
4. **Calibrate the receivers.** Send 64 known input vectors and fit each output channel's LO phase and gain (digital).
5. **Accept** if every core reaches ENOB ≥ 4 against held-out vectors [simulated: 4.8].

## 7. Run

* The host streams input vectors to the modulators [64 channels at 10 GBd] and reads the coherent receivers. Each symbol performs one 16×16 matrix–vector product per core [20.5 TOPS].
* New weights: re-expose only the tiles that contain weight patches (listed in `tiles.json`, `contains_weights`). For stationary-weight inference this is rare. For training-style updates it is the bottleneck (ms–s per tile).
* Monitor drift: re-run step 6.4 (receiver phases) every few minutes and step 6.1–6.3 when the fidelity on reference vectors drops.

## 8. Safety

* 532 nm writing lasers at ≥ 1 W are Class 4: enclose the beam path, interlock the enclosure, and wear OD ≥ 6 eyewear at 532 nm during alignment.
* 1550 nm readout lasers above 10 mW are Class 3B: never view fiber ends.
* Sb₂Se₃ contains selenium: handle targets and sputter-chamber waste as toxic.

## 9. What is still open

* Long-term stability of multi-level Sb₂Se₃ states, and their drift with temperature.
* Thermal crosstalk between neighbouring weight patches while writing (patch pitch ≈ 13 µm).
* The vertical mode and facet coupling of a 20 nm film on 300 nm SiN: verify with a 3-D mode solver.
* Wavelength dependence of the couplers. This sets how many WDM channels can share a core.
* Coupler spread limits precision to ≈5 bits. Redundant MZI designs (Hamerly et al. 2022) are the known route to more.
