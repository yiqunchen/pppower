// Run the browser calculator in vignettes/sample-size-calculator.Rmd under
// Node with a minimal DOM stub, and print the displayed results as JSON.
// Usage: node calculator_harness.js '{"dm_r2_A": 0.3}'
const fs = require("fs");
const path = require("path");
const rmd = fs.readFileSync(path.join(__dirname, "../../vignettes/sample-size-calculator.Rmd"), "utf8");
const overrides = JSON.parse(process.argv[2] || "{}");

const els = {};
function el(idName) {
  if (!els[idName]) {
    els[idName] = {
      id: idName, value: "", textContent: "", style: {}, handlers: {},
      classList: { add() {}, remove() {}, toggle() {} },
      addEventListener(ev, fn) { (this.handlers[ev] = this.handlers[ev] || []).push(fn); },
      setAttribute() {}, getAttribute() { return null; }
    };
  }
  return els[idName];
}
for (const m of rmd.matchAll(/<input id="([^"]+)"[^>]*value="([^"]*)"/g)) el(m[1]).value = m[2];
for (const m of rmd.matchAll(/<select id="([^"]+)">([\s\S]*?)<\/select>/g)) {
  const sel = m[2].match(/<option value="([^"]+)" selected/);
  el(m[1]).value = sel ? sel[1] : (m[2].match(/<option value="([^"]+)"/) || [])[1];
}
for (const m of rmd.matchAll(/id="([^"]+)"/g)) el(m[1]);

global.document = {
  readyState: "complete",
  getElementById: (x) => els[x] || null,
  querySelectorAll: () => [],
  addEventListener() {}
};
const script = rmd.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
eval(script);

for (const [k, v] of Object.entries(overrides)) {
  el(k).value = String(v);
  (el(k).handlers.input || []).forEach((fn) => fn());
}
const out = {};
for (const p of ["single", "diff", "dm", "pm", "pb"]) {
  out[p] = {};
  for (const suffix of ["_n_ppi", "_n_classic", "_n_ppi_group", "_n_classic_group", "_message"]) {
    if (els[p + suffix]) out[p][suffix.slice(1)] = els[p + suffix].textContent;
  }
}
console.log(JSON.stringify(out));
