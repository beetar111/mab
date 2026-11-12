#!/usr/bin/env node
// Inline core/otc.js and the default design into designer/index.html (single self-contained file).
const fs = require("fs"), path = require("path");
const root = path.join(__dirname, "..");
const design = process.argv[2] || path.join(root, "designs/ogpu-1.json");
let html = fs.readFileSync(path.join(root, "designer/src.html"), "utf8");
html = html.replace("/*__OTC__*/", () => fs.readFileSync(path.join(root, "core/otc.js"), "utf8"));
html = html.replace("/*__DESIGN__*/null", () => JSON.stringify(JSON.parse(fs.readFileSync(design, "utf8"))));
fs.writeFileSync(path.join(root, "designer/index.html"), html);
console.log("designer/index.html", (html.length / 1024).toFixed(0), "KB");
