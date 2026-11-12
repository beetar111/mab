#!/usr/bin/env node
/* Compile an optical-GPU design file into manufacturing outputs.
 *   node tools/run.js designs/ogpu-1.json build/ogpu-1 [--all-tiles] [--holo 2] [--realizations 10]
 * Outputs: resolved design + weight program, GDSII layout, DMD frames (programmable plate), LCoS phase holograms
 * (Fresnel plate), overview images, report.md.
 */
const fs = require("fs"), path = require("path"), zlib = require("zlib");
const OTC = require("../core/otc.js");

const args = process.argv.slice(2);
const designPath = args[0] || path.join(__dirname, "../designs/ogpu-1.json");
const outDir = args[1] || path.join(__dirname, "../build", path.basename(designPath, ".json"));
const flag = (n, d) => { const i = args.indexOf(n); return i < 0 ? d : (args[i + 1] && !args[i + 1].startsWith("--") ? args[i + 1] : true); };
const ALL = !!flag("--all-tiles", false), HOLO = +flag("--holo", 2), REAL = +flag("--realizations", 10);
const deflate = async u => new Uint8Array(zlib.deflateSync(Buffer.from(u), { level: 9 }));
const mk = d => fs.mkdirSync(d, { recursive: true });
const fmt = (v, d = 2) => Number(v).toFixed(d);
const t0 = Date.now();

(async () => {
  mk(outDir); mk(path.join(outDir, "plate/dmd")); mk(path.join(outDir, "plate/lcos")); mk(path.join(outDir, "program"));
  const design = JSON.parse(fs.readFileSync(designPath, "utf8"));
  const R = await OTC.run(design, { realizations: REAL, samples: 300, log: (m) => process.stdout.write(`  ${m}\n`) });
  const { g, layout, perf, acc, progs, proc, stats } = R;

  // ── weight program (what the writer exposes as graded patches) ──
  const plat = OTC.PLATFORMS[design.platform];
  for (let ci = 0; ci < design.cores.length; ci++) {
    const rows = ["block,col,mode,kind,x0_um,x1_um,y_um,dose_level_0_255,phase_rad"];
    for (const p of layout.patches.filter(q => q.core === ci)) rows.push([p.block, p.col, p.m, p.kind, fmt(p.x0, 3), fmt(p.x1, 3), fmt(p.y, 3), Math.round(255 * p.level), fmt(p.value, 5)].join(","));
    fs.writeFileSync(path.join(outDir, `program/core${ci}_weights.csv`), rows.join("\n") + "\n");
  }
  // ── resolved design ──
  const resolved = {
    design, generated: new Date().toISOString(), engine: "OTC " + OTC.version,
    cell: { core_width_um: g.w, lane_pitch_um: g.D, coupler_gap_um: g.gap, coupler_length_um: g.L50, sbend_length_um: g.Ls, phase_patch_length_um: g.Lp, mzi_length_um: g.Lcell, min_bend_radius_um: g.Rmin, neff: g.neff, write_blur_sigma_um: g.cfg.blur, speckle_rms: g.cfg.sigF },
    process_centering: g.centering,
    device_pool: stats.pool,
    die_um: layout.die, tiles: { plate: proc.plate, ...R.tiles, slm: undefined },
    accuracy: acc.map((a, i) => ({ core: i, weights: design.cores[i].weights, fidelity_calibrated: a.cal.fid, enob_calibrated: a.cal.run.enob, fidelity_uncalibrated: a.uncal.fid, enob_uncalibrated: a.uncal.run.enob })),
    performance: perf,
  };
  fs.writeFileSync(path.join(outDir, "design.resolved.json"), JSON.stringify(resolved, (k, v) => (typeof v === "number" ? +v.toPrecision(6) : v), 2));

  // ── GDSII ──
  const gds = OTC.gds(layout, design.name.replace(/[^A-Za-z0-9_]/g, "_"));
  fs.writeFileSync(path.join(outDir, "layout.gds"), gds);
  console.log(`  GDSII: ${(gds.length / 1e6).toFixed(2)} MB`);

  // ── overview image (whole die, 10 µm pixels) and a zoom on the first MZI columns ──
  const ovPx = Math.max(5, layout.die.W / 8000);
  const ovW = Math.ceil(layout.die.W / ovPx), ovH = Math.ceil(layout.die.H / ovPx);
  const thick = Object.assign({}, layout, { paths: layout.paths.map(p => Object.assign({}, p, { w: Math.max(p.w, ovPx * 1.2) })), patches: layout.patches.map(q => Object.assign({}, q, { w: Math.max(q.w, ovPx * 1.2) })) });
  fs.writeFileSync(path.join(outDir, "overview.png"), await OTC.pngGray(OTC.raster(thick, 0, 0, ovW, ovH, ovPx), ovW, ovH, deflate));
  const blk = layout.blocks.find(b => b.core === 0 && b.type === "meshV");
  const zx = blk.x0 + g.Lphi - 50, zy = blk.y0 + (blk.y1 - blk.y0) / 2 - 70, zp = g.w / 4;
  const zW = Math.round(Math.min(g.Lcell, 1000) / zp), zH = Math.round(140 / zp);
  fs.writeFileSync(path.join(outDir, "zoom_mzi.png"), await OTC.pngGray(OTC.raster(layout, zx, zy, zW, zH, zp), zW, zH, deflate));

  // ── programmable plate A: DMD frames (direct imaging, 8-bit gray by PWM) ──
  const dmdPlan = OTC.tilePlan(layout, Object.assign({}, proc, { plate: "dmd" }));
  const [dw, dh] = proc.dmd.pixels, tiles = [];
  let written = 0;
  for (let r = 0; r < dmdPlan.ny; r++) for (let c = 0; c < dmdPlan.nx; c++) {
    const x0 = c * dmdPlan.tw, y0 = r * dmdPlan.th;
    const hasW = layout.patches.some(q => q.x1 > x0 && q.x0 < x0 + dmdPlan.tw && q.y > y0 && q.y < y0 + dmdPlan.th);
    const hasG = layout.paths.some(p => p.pts.some(q => q[0] >= x0 && q[0] < x0 + dmdPlan.tw && q[1] >= y0 && q[1] < y0 + dmdPlan.th));
    tiles.push({ row: r, col: c, x0_um: +x0.toFixed(2), y0_um: +y0.toFixed(2), contains_waveguides: hasG, contains_weights: hasW });
    const sample = (r === Math.floor(dmdPlan.ny / 2) && c < 3) || (r === 0 && c === 0);
    if ((ALL && hasG) || sample) {
      const img = OTC.raster(layout, x0, y0, dw, dh, dmdPlan.px);
      fs.writeFileSync(path.join(outDir, `plate/dmd/tile_r${String(r).padStart(2, "0")}_c${String(c).padStart(3, "0")}.png`), await OTC.pngGray(img, dw, dh, deflate));
      written++;
    }
  }
  fs.writeFileSync(path.join(outDir, "plate/dmd/tiles.json"), JSON.stringify({ device: proc.dmd, slab_pixel_um: dmdPlan.px, demagnification: dmdPlan.demag, tile_um: [dmdPlan.tw, dmdPlan.th], grid: [dmdPlan.nx, dmdPlan.ny], tiles }, null, 1));
  console.log(`  DMD: ${dmdPlan.nx}×${dmdPlan.ny} tiles planned, ${written} frames written${ALL ? "" : " (samples; use --all-tiles for the full set)"}`);

  // ── programmable plate B / static Fresnel plate: LCoS phase holograms (GS, 1024² sub-tiles) ──
  const HN = 1024, hp = g.w / 4, holo = [];
  const centers = [[blk.x0 + 0.5 * g.Lcell, (blk.y0 + blk.y1) / 2], [layout.blocks.find(b => b.core === 0 && b.type === "sigma").x0 + 0.6 * g.Lcell, (blk.y0 + blk.y1) / 2]];
  for (let k = 0; k < Math.min(HOLO, centers.length); k++) {
    const fill = 0.7, n = Math.round(HN * fill), x0 = centers[k][0] - n * hp / 2, y0 = centers[k][1] - n * hp / 2;
    const sub = OTC.raster(layout, x0, y0, n, n, hp), T = new Float64Array(HN * HN), o = Math.floor((HN - n) / 2);
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) T[(j + o) * HN + i + o] = sub[j * n + i] / 255;
    const t1 = Date.now(), H = OTC.hologram(T, HN, { iters: 20, fill: 0.75, seed: 3 + k });
    const plate = OTC.fresnelPlate(H.phase, HN, 256);
    let mx = 0; for (let p = 0; p < HN * HN; p++) mx = Math.max(mx, H.I[p]);
    const rec = new Uint8Array(HN * HN); for (let p = 0; p < HN * HN; p++) rec[p] = Math.min(255, Math.round(255 * Math.sqrt(H.I[p] / (H.scale * 1.3))));
    const tgt = new Uint8Array(HN * HN); for (let p = 0; p < HN * HN; p++) tgt[p] = Math.round(255 * T[p]);
    fs.writeFileSync(path.join(outDir, `plate/lcos/holo${k}_phase_fresnel.png`), await OTC.pngGray(plate, HN, HN, deflate));
    fs.writeFileSync(path.join(outDir, `plate/lcos/holo${k}_target.png`), await OTC.pngGray(tgt, HN, HN, deflate));
    fs.writeFileSync(path.join(outDir, `plate/lcos/holo${k}_reconstruction.png`), await OTC.pngGray(rec, HN, HN, deflate));
    holo.push({ k, x0_um: x0, y0_um: y0, size_px: HN, slab_pixel_um: hp, speckle_rms_single_frame: H.sigma, seconds: (Date.now() - t1) / 1000 });
    console.log(`  LCoS hologram ${k}: speckle ${fmt(100 * H.sigma, 1)}% rms (single frame), ${fmt((Date.now() - t1) / 1000, 1)} s`);
  }
  const f_mm = HN * proc.lcos.pitch_um ** 2 / (plat.writeNm / 1000) / 1000;
  fs.writeFileSync(path.join(outDir, "plate/lcos/holograms.json"), JSON.stringify({ size: HN, plate_pitch_um: proc.lcos.pitch_um, lensless_distance_mm: f_mm, note: "phase 0..255 ↔ 0..2π, includes Fresnel lens term; display on a phase-only SLM or etch as a 256-level DOE", holo }, null, 1));

  // ── report ──
  const P = perf, w = s => s.toFixed(2);
  const tilesLcos = Math.ceil(layout.die.W / (HN * hp)) * Math.ceil(layout.die.H / (HN * hp));
  const writeTime = n => n * (proc.exposure_s + proc.settle_s) / 60;
  const report = `# ${design.name} — compiled design report

Generated by \`tools/run.js\` (engine OTC ${OTC.version}) from \`${path.relative(process.cwd(), designPath)}\` in ${fmt((Date.now() - t0) / 1000, 0)} s.

${design.description}

## Platform and cell
| | |
|---|---|
| Platform | ${plat.name} |
| Layer stack | ${plat.stack.join(" → ")} |
| Index contrast used for guides / weights | ${(plat.dn * plat.wgLevel).toExponential(1)} / ${(plat.dn * (1 - plat.wgLevel)).toExponential(1)} |
| Core width, lane pitch | ${fmt(g.w, 2)} µm, ${fmt(g.D, 1)} µm |
| Directional coupler (gap, length) | ${fmt(g.gap, 2)} µm, ${fmt(g.L50, 0)} µm (process-centered: mean α ${g.centering.map(c => fmt(c.alphaMean, 4)).join(" → ")} rad, target 0.7854) |
| Weight patch length (2π) | ${fmt(g.Lp, 0)} µm |
| MZI cell length | ${fmt(g.Lcell, 0)} µm |
| Minimum bend radius (0.1 dB/cm) | ${fmt(g.Rmin / 1000, 2)} mm |
| Written-index blur σ, speckle | ${fmt(g.cfg.blur, 2)} µm, ${fmt(100 * g.cfg.sigF, 1)} % rms (${proc.plate.toUpperCase()}, ${proc.frames} frames) |
| Measured MZI pool (${stats.pool.length} BPM realizations) | α₁ ${fmt(Math.min(...stats.pool.map(p => p.a1)), 3)}–${fmt(Math.max(...stats.pool.map(p => p.a1)), 3)} rad, IL ${fmt(stats.pool.reduce((a, p) => a + p.IL, 0) / stats.pool.length, 2)} dB, cross extinction ${fmt(Math.min(...stats.pool.map(p => p.erCross)), 1)}–${fmt(Math.max(...stats.pool.map(p => p.erCross)), 1)} dB |

## Accuracy (each core: ${design.N}×${design.N}, 300 random input vectors, ${design.system.dacBits}-bit DAC, ${design.system.adcBits}-bit ADC, receiver SNR ${fmt(R.rx.snrDb, 1)} dB)
| Core | Weights | Fidelity, calibrated | ENOB, calibrated | Fidelity, uncalibrated | ENOB, uncalibrated |
|---|---|---|---|---|---|
${acc.map((a, i) => `| ${i} | ${design.cores[i].weights} | ${fmt(a.cal.fid, 4)} | ${fmt(a.cal.run.enob, 2)} | ${fmt(a.uncal.fid, 4)} | ${fmt(a.uncal.run.enob, 2)} |`).join("\n")}

"Calibrated" = every MZI measured after writing, phases re-solved for the real couplers (hardware-aware programming) and written closed-loop to 1/${plat.levels} of 2π. Uncalibrated meshes do not compute; calibration is part of the process, not an option.

## Performance (model; assumptions in the design file)
| | This chip (light-written weights) | Same chip, silicon-photonics thermo-optic weights | ${P.gpu.name} |
|---|---|---|---|
| Throughput | ${fmt(P.ops / 1e12, 1)} TOPS (${P.K} cores × ${P.N}² MAC × ${design.system.symbolRate_GHz} GHz${P.lambdas > 1 ? ` × ${P.lambdas} λ` : ""}) | same | ${P.gpu.tops} TOPS |
| Power | ${w(P.total)} W | ${w(P.si.total)} W (weights hold ${w(P.si.weights)} W) | ${P.gpu.watts} W |
| Efficiency | **${fmt(P.topsW, 1)} TOPS/W** | ${fmt(P.si.topsW, 2)} TOPS/W | ${fmt(P.gpu.topsW, 2)} TOPS/W |
| Precision | ≈${fmt(P.enob, 1)} effective bits | similar | 8-bit integer exact |
| Latency through the die | ${fmt(P.latency_ns, 2)} ns | | |

Power breakdown: laser ${w(P.P.laser)} W, modulators ${w(P.P.modulators)} W, DACs ${w(P.P.dacs)} W, receivers + ADCs ${w(P.P.receivers)} W, weight hold ${w(P.P.weights)} W, temperature control ${w(P.P.thermal)} W.
The comparison is per useful operation at different precision; the optical core is a ~5-bit analog matrix engine attached to an electronic host, not a drop-in GPU.

## Die and writing
| | |
|---|---|
| Die | ${fmt(layout.die.W / 1000, 1)} mm × ${fmt(layout.die.H / 1000, 1)} mm (${P.nMzi} MZIs incl. Σ, ${layout.patches.length} weight patches) |
| Slab pixel (w/4) | ${fmt(g.w / 4, 3)} µm |
| DMD direct-imaging tiles | ${dmdPlan.nx} × ${dmdPlan.ny} = ${dmdPlan.count} (demag ×${fmt(dmdPlan.demag, 1)}), full write ≈ ${fmt(writeTime(tiles.filter(t => t.contains_waveguides).length), 1)} min |
| LCoS holographic tiles (1024² each) | ${tilesLcos}, Fresnel plate distance ${fmt(f_mm, 1)} mm |
| Weight-only rewrite | ${tiles.filter(t => t.contains_weights).length} DMD tiles ≈ ${fmt(writeTime(tiles.filter(t => t.contains_weights).length), 1)} min |

## Files
- \`design.resolved.json\`: design + derived geometry, device statistics, accuracy, performance
- \`layout.gds\`: GDSII (1 nm grid). Layer 1 waveguides (datatype = write dose 0–255), 2 weight patches (datatype = dose), 3 facet converters / terminators, 10 fiducials, 63 die outline
- \`program/core*_weights.csv\`: weight patch list for the write–verify loop
- \`plate/dmd/\`: DMD frames (8-bit, one per exposure tile) + \`tiles.json\` stage map
- \`plate/lcos/\`: Fresnel phase plates (256-level), targets and simulated reconstructions
- \`overview.png\`, \`zoom_mzi.png\`: raster previews

See \`docs/BUILD.md\` for how to build and operate it.
`;
  fs.writeFileSync(path.join(outDir, "report.md"), report);
  console.log(`Done in ${fmt((Date.now() - t0) / 1000, 0)} s → ${outDir}`);
})().catch(e => { console.error(e); process.exit(1); });
