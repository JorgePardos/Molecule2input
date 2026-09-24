/* m2i web interface.
 *
 * Plain JavaScript, no build step. The page is re-rendered from `state` after
 * every committed change (a file, a SMILES, a choice); typing in a field only
 * updates `state`, so the caret is never lost. Everything a chemist reads as
 * data comes from the API as it is -- the engine's own words, shown verbatim.
 * Text from the server or from files is always set as text, never as HTML.
 */
"use strict";

(() => {
  // -- tiny DOM helper ----------------------------------------------------------------
  function h(tag, attrs, ...children) {
    const el = document.createElement(tag);
    for (const [key, value] of Object.entries(attrs || {})) {
      if (value === null || value === undefined || value === false) continue;
      if (key === "class") el.className = value;
      else if (key === "style" && typeof value === "object") Object.assign(el.style, value);
      else if (key.startsWith("on") && typeof value === "function") el.addEventListener(key.slice(2), value);
      else if (key === "dataset") Object.assign(el.dataset, value);
      else if (value === true) el.setAttribute(key, "");
      else if (key in el && !["list", "form", "type"].includes(key)) el[key] = value;
      else el.setAttribute(key, value);
    }
    for (const child of children.flat(Infinity)) {
      if (child === null || child === undefined || child === false) continue;
      el.append(child instanceof Node ? child : document.createTextNode(String(child)));
    }
    return el;
  }
  const $ = (id) => document.getElementById(id);

  async function api(method, url, body, form) {
    const options = { method, credentials: "same-origin", headers: {} };
    if (form) options.body = form;
    else if (body !== undefined) {
      options.body = JSON.stringify(body);
      options.headers["Content-Type"] = "application/json";
    }
    let response;
    try {
      response = await fetch(url, options);
    } catch (err) {
      throw new Error("The server could not be reached. Is it still running?");
    }
    const data = await response.json().catch(() => ({}));
    if (!response.ok) throw new Error(data.error || `The server answered ${response.status}.`);
    return data;
  }

  // -- state --------------------------------------------------------------------------------
  const KINDS = [
    { id: "picture", name: "Picture", body: "Read by DECIMER. Confirmed when in doubt." },
    { id: "file", name: "ChemDraw / molfile", body: "Exactly, wedges included. Nothing asked." },
    { id: "smiles", name: "SMILES", body: "Typed. Exactly as written." },
    { id: "cif", name: "Crystal (.cif)", body: "Rebuilt with its measured geometry." },
  ];
  const FORMATS = { gaussian: [".gjf", "Gaussian"], orca: [".inp", "ORCA"], xyz: [".xyz", "XYZ"], sdf: [".sdf", "SDF"] };
  const FORMAT_LABEL = { gaussian: "Gaussian (.gjf)", orca: "ORCA (.inp)", xyz: "XYZ (.xyz)", sdf: "SDF (.sdf)" };

  const state = {
    about: null,
    health: null,
    splashDismissed: false,
    kind: "picture",
    sources: { picture: null, file: null, smiles: null, cif: null },
    checks: {},
    errors: {},
    busy: {},
    expanded: false,
    notThis: false,
    advanced: { charge: "", mult: "", keepAll: false, keep: 1, nconfs: "", forceField: "mmff94s", seed: 61453 },
    output: { format: "gaussian", recipes: {}, recipe: {}, fields: {}, problem: null, jobs: "", name: "" },
    result: null,
    smilesDraft: "",
    om: null,
    cif: null,
  };

  // -- boot ------------------------------------------------------------------------------------
  async function boot() {
    try {
      state.about = await api("GET", "/api/about");
    } catch (err) {
      state.errors.boot = err.message;
    }
    await pollHealth();
    await loadRecipes("gaussian");
    render();
  }

  async function pollHealth() {
    try {
      state.health = await api("GET", "/api/health");
    } catch (err) {
      state.health = null;
    }
    const model = state.health && state.health.model.state;
    if (model === "waiting" || model === "loading") {
      setTimeout(async () => { await pollHealth(); renderSplash(); renderModelLine(); }, 2000);
    } else if (!$("splash").hidden) {
      state.splashDismissed = true; // the page comes in on its own
      render();
    }
  }

  // -- rendering ---------------------------------------------------------------------------------
  function render() {
    const focus = rememberFocus();
    renderHeader();
    renderNav();
    const main = $("main");
    main.replaceChildren(...renderMain().flat(Infinity).filter(Boolean));
    renderSplash();
    restoreFocus(focus);
  }

  function rememberFocus() {
    const el = document.activeElement;
    if (!el || !el.dataset || !el.dataset.focus) return null;
    return { id: el.dataset.focus, start: el.selectionStart, end: el.selectionEnd };
  }

  function restoreFocus(focus) {
    if (!focus) return;
    const el = document.querySelector(`[data-focus="${CSS.escape(focus.id)}"]`);
    if (!el) return;
    el.focus();
    if (focus.start !== undefined && focus.start !== null && el.setSelectionRange) {
      try { el.setSelectionRange(focus.start, focus.end); } catch (err) { /* not a text field */ }
    }
  }

  function currentSource() {
    return state.sources[state.kind];
  }

  function currentCheck() {
    const source = currentSource();
    return source && source.key ? state.checks[source.key] : null;
  }

  function renderHeader() {
    const source = currentSource();
    let text;
    if (state.busy.reading) text = [h("span", { class: "spinner small" }), "reading"];
    else if (state.result && resultCurrent()) text = "ready";
    else if (state.kind === "file" && state.om) text = "metal complex · from a drawing";
    else if (state.kind === "cif" && state.cif && state.cif.summary) text = "crystal · rebuilt as measured";
    else if (source && state.kind === "picture" && source.read) text = "picture · decimer";
    else if (source && state.kind === "smiles") text = "smiles";
    else if (source && state.kind === "file") text = "chemdraw";
    else {
      const version = state.about ? `v${state.about.version}` : "";
      text = state.about && state.about.hosted
        ? [version, h("span", { class: "long" }, " · private session · files removed a day after last use")]
        : [version, h("span", { class: "long" }, " · running on this computer")];
    }
    $("context").replaceChildren(...[text].flat());
  }

  // -- navigation spine ------------------------------------------------------------------------------
  function steps() {
    const source = currentSource();
    const check = currentCheck();
    const result = state.result && resultCurrent();
    if (state.kind === "file" && state.om) return complexSteps();
    if (state.kind === "cif" && state.cif && state.cif.summary) return crystalSteps();

    const given = source && (source.read || state.kind !== "picture");
    const structureSub = !source ? "nothing given yet"
      : state.busy.reading ? "reading the photo…"
      : state.kind === "smiles" ? "smiles · typed"
      : state.kind === "file" ? `${source.shown_name || "file"}`
      : source.read ? "picture · decimer" : "photo given";
    const s1 = { name: "Structure", sub: structureSub, status: given ? "done" : "active", subClass: state.busy.reading ? "busy" : "" };

    let s2 = { name: "Check", status: "pending" };
    let passed = false;
    if (given && check) {
      if (check.molecule === null || check.errors) {
        s2 = { name: "Check", sub: "needs correcting", subClass: "warn", status: "active" };
      } else if (!check.gate.needed) {
        s2 = { name: "Check", sub: "not checked", subClass: "good", status: "done" };
        passed = true;
      } else if (check.gate.confirmed) {
        s2 = { name: "Check", sub: "confirmed by you", subClass: "good", status: "done" };
        passed = true;
      } else {
        const n = check.gate.reasons.length;
        s2 = { name: "Check", sub: `${n} reason${n === 1 ? "" : "s"} to look`, subClass: "warn", status: "active" };
      }
    }
    const s3 = { name: "Output", status: result ? "done" : passed ? "active" : "pending",
      sub: result ? `${state.result.format} · ${state.result.jobs}` : passed ? "choose a format" : null };
    return [s1, s2, s3];
  }

  function renderNav() {
    const list = steps();
    const nav = $("nav");
    nav.replaceChildren(...[
      ...list.map((step, i) => h("div", { class: `step ${step.status !== "pending" ? "on" : ""}` },
        h("div", { class: "step-track" },
          h("span", { class: `step-dot ${step.status}` }, step.status === "done" ? "✓" : String(i + 1)),
          i < list.length - 1 ? h("span", { class: "step-line" }) : null),
        h("div", null,
          h("div", { class: "step-name" }, step.name),
          step.sub ? h("div", { class: `step-sub ${step.subClass || ""}` }, step.sub) : null))),
      showAdvanced() ? renderAdvanced() : null,
    ].filter(Boolean));
    $("progress").replaceChildren(...list.map((step) => h("div", { class: step.status })));
  }

  function showAdvanced() {
    return !(state.kind === "cif" || (state.kind === "file" && state.om));
  }

  function renderAdvanced() {
    const a = state.advanced;
    const set = (key, recheck) => (event) => {
      const input = event.target;
      a[key] = input.type === "checkbox" ? input.checked : input.value;
      if (recheck) refreshCheck();
      else render();
    };
    const summary = (v, fallback) => (v === "" || v === null ? fallback : v);
    return h("details", { class: "advanced", open: state.advancedOpen, ontoggle: (e) => { state.advancedOpen = e.target.open; } },
      h("summary", null, "Advanced"),
      state.advancedOpen ? [
        advRow("Charge", h("input", { type: "number", step: 1, placeholder: "auto", value: a.charge, "data-focus": "adv-charge", onchange: set("charge", true) })),
        advRow("Multiplicity", h("input", { type: "number", min: 1, step: 1, placeholder: "auto", value: a.mult, "data-focus": "adv-mult", onchange: set("mult", true) })),
        advRow("Conformers kept", h("input", { type: "number", min: 1, max: 10, step: 1, value: a.keep, onchange: set("keep") })),
        advRow("Embedded", h("input", { type: "number", min: 0, step: 10, placeholder: "auto", value: a.nconfs, onchange: set("nconfs") })),
        advRow("Force field", h("select", { onchange: set("forceField") },
          ["mmff94s", "mmff94", "uff", "none"].map((ff) => h("option", { value: ff, selected: ff === a.forceField }, ff)))),
        advRow("Random seed", h("input", { type: "number", step: 1, value: a.seed, onchange: set("seed") })),
        h("div", { class: "adv-row" }, h("label", { class: "check" },
          h("input", { type: "checkbox", checked: a.keepAll, onchange: set("keepAll", true) }), "Keep every fragment")),
      ] : [
        advSummary("Charge", summary(a.charge, "auto")),
        advSummary("Multiplicity", summary(a.mult, "auto")),
        advSummary("Conformers", `${a.keep} / ${summary(a.nconfs, "auto")}`),
        advSummary("Force field", a.forceField),
      ]);
  }
  const advRow = (label, input) => h("div", { class: "adv-row" }, h("span", null, label), input);
  const advSummary = (label, value) => h("div", { class: "adv-row" }, h("span", null, label), h("span", { class: "value" }, String(value)));

  // -- main --------------------------------------------------------------------------------------------
  function renderMain() {
    if (state.errors.boot) return [strip("error", state.errors.boot)];
    const out = [renderStructureStep()];
    if (state.kind === "file" && state.om) out.push(...renderComplex());
    else if (state.kind === "cif") out.push(...renderCrystal());
    else out.push(...renderOrganic());
    return out;
  }

  function nothingGiven() {
    return !Object.values(state.sources).some(Boolean) && !state.om && !(state.cif && state.cif.summary);
  }

  function renderStructureStep() {
    const section = h("section", { class: "section" }, h("div", { class: "eyebrow" }, "Step 1 · Structure"));
    if (nothingGiven()) {
      section.append(
        h("h2", { class: "title" }, "Start from"),
        h("p", { class: "intro" }, "A photo of a hand-drawn molecule, a ChemDraw file or molfile, a SMILES, or a crystal structure. If you have the ChemDraw file, use it rather than a picture of it: the file holds the bonds, wedges and labels exactly."),
        h("div", { class: "cards", role: "radiogroup", "aria-label": "Start from" },
          KINDS.map((kind) => h("button", {
            class: `card ${state.kind === kind.id ? "selected" : ""}`, role: "radio", "aria-checked": String(state.kind === kind.id),
            onclick: () => chooseKind(kind.id),
          }, h("div", { class: "card-title" }, kind.name), h("div", { class: "card-body" }, kind.body)))),
      );
    } else {
      section.append(h("div", { class: "segmented", role: "tablist" },
        KINDS.map((kind) => h("button", { class: state.kind === kind.id ? "on" : "", role: "tab", "aria-selected": String(state.kind === kind.id), onclick: () => chooseKind(kind.id) }, kind.name))));
    }
    section.append(...renderSourceInput().filter(Boolean));
    return section;
  }

  function chooseKind(kind) {
    state.kind = kind;
    state.errors.source = null;
    state.notThis = false;
    render();
  }

  function renderSourceInput() {
    const source = currentSource();
    const types = (state.about && state.about.types) || {};
    const limit = state.about ? state.about.max_upload_mb : 20;
    const error = state.errors.source ? strip("error", state.errors.source) : null;
    if (state.kind === "smiles") {
      return [
        h("input", {
          class: "smiles-input", placeholder: "C[C@H](N)C(=O)O", spellcheck: false, autocomplete: "off", "aria-label": "SMILES",
          value: state.smilesDraft, "data-focus": "smiles",
          oninput: (e) => { state.smilesDraft = e.target.value; },
          onkeydown: (e) => { if (e.key === "Enter") { state.smilesDraft = e.target.value; submitSmiles(); } },
          onchange: (e) => { state.smilesDraft = e.target.value; submitSmiles(); },
        }),
        error,
      ];
    }
    if (state.kind === "picture") {
      if (!source) {
        return [dropzone("Drop a photo here, or", types.picture || [], limit, uploadPicture, "image/*"), error, modelLine()];
      }
      if (state.busy.reading) return [renderReading(source), error];
      const row = fileRow(source.preview_url ? h("img", { src: source.preview_url, alt: "" }) : "IMG",
        source.shown_name, source.read ? `read by decimer · ${source.style === "hand_drawn" ? "hand-drawn" : "clean drawing"}` : "not read yet",
        () => { state.sources.picture = null; state.result = null; render(); });
      const again = !source.read ? h("div", { class: "row", style: { marginTop: "10px" } },
        h("button", { class: "btn primary", onclick: () => readPicture(source) }, "Read the structure")) : null;
      return [row, error, again, source.read ? null : modelLine()];
    }
    if (state.kind === "file") {
      if (!source && !state.om) {
        return [dropzone("Drop a ChemDraw file or molfile here, or", types.file || [], limit, uploadFile, ".cdx,.cdxml,.mol,.sdf,.mdl"), error];
      }
      const shown = state.om ? state.om.shown_name : source.shown_name;
      const meta = state.om ? "metal complex · bonds to the metal read as drawn"
        : source.atoms !== undefined ? `${source.atoms} atoms · ${source.bonds} bonds · ${source.wedges} wedge${source.wedges === 1 ? "" : "s"}` : "read from the file";
      return [
        fileRow("CDX", shown, meta, () => { state.sources.file = null; state.om = null; state.result = null; render(); }),
        h("div", { class: "caption", style: { marginBottom: "4px" } }, "Read straight from the file: bonds and wedges are exactly as drawn, with no recognition model involved."),
        error,
      ];
    }
    // crystal
    if (!state.cif || !state.cif.summary) {
      return [dropzone("Drop a crystal structure (.cif) here, or", ["cif"], limit, uploadCif, ".cif"), error,
        h("div", { class: "caption", style: { marginTop: "12px" } }, "The geometry is taken from the crystal as measured, so stereochemistry and ligand shapes are kept exactly. Oxidation state, charge and spin are not in a CIF: you will be asked for them.")];
    }
    return [renderCrystalFile(), error];
  }

  function dropzone(prompt, types, limit, onfile, accept) {
    const input = h("input", { type: "file", accept, onchange: (e) => { if (e.target.files[0]) onfile(e.target.files[0]); } });
    const zone = h("label", {
      class: "dropzone",
      ondragover: (e) => { e.preventDefault(); zone.classList.add("over"); },
      ondragleave: () => zone.classList.remove("over"),
      ondrop: (e) => { e.preventDefault(); zone.classList.remove("over"); if (e.dataTransfer.files[0]) onfile(e.dataTransfer.files[0]); },
    },
    h("img", { src: "logo.svg", alt: "" }),
    h("div", { class: "dropzone-title" }, `${prompt} `, h("u", null, "browse")),
    h("div", { class: "dropzone-types" }, `${types.join(" · ")} · up to ${limit} MB`),
    input);
    if (state.busy.upload) return h("div", { class: "dropzone" }, h("div", { class: "row", style: { justifyContent: "center" } }, h("span", { class: "spinner" }), "Uploading…"));
    return zone;
  }

  function fileRow(badge, name, meta, onreplace) {
    return h("div", { class: "file-row" },
      h("div", { class: "file-badge" }, badge),
      h("div", { class: "grow" }, h("div", { class: "file-name" }, name || ""), h("div", { class: "file-meta" }, meta)),
      h("button", { class: "btn ghost", onclick: onreplace }, "Replace"));
  }

  function modelLine() {
    const model = state.health && state.health.model;
    let dot = "grey", text = "checking the recognition model…";
    if (model) {
      if (model.state === "ready") { dot = ""; text = `${model.backend} loaded and ready — best at hand-drawn structures`; }
      else if (model.state === "loading" || model.state === "waiting") { dot = "amber"; text = `loading the recognition model (${model.elapsed} s of about ${model.typical} s) — a photo can be given meanwhile`; }
      else if (model.state === "failed") { dot = "red"; text = `the recognition model could not start: ${model.error}`; }
      else { dot = "red"; text = "no recognition model on this server: type the SMILES, or upload the ChemDraw file"; }
    }
    return h("div", { class: "model-line", id: "model-line" }, h("span", { class: `dot ${dot}` }), text);
  }

  function renderModelLine() {
    const line = $("model-line");
    if (line) line.replaceWith(modelLine());
  }

  function renderReading(source) {
    const model = state.health && state.health.model;
    const loaded = model && model.state === "ready";
    return h("div", { class: "reading" },
      h("div", { class: "panel pane" },
        h("div", { class: "pane-title" }, source.shown_name || "photo"),
        h("div", { class: "pane-box" }, source.preview_url ? h("img", { src: source.preview_url, alt: "The photo being read" }) : null)),
      h("div", { class: "reading-panel", role: "status", "aria-live": "polite" },
        h("div", { class: "row" }, h("span", { class: "spinner" }), h("div", { style: { font: "600 15px var(--sans)" } }, "Reading the drawing")),
        h("div", { class: "caption", style: { marginTop: "10px" } }, loaded
          ? "DECIMER is already in memory, so this is a couple of seconds. Only the first reading after a start is the slow one."
          : "The recognition model is still loading; the reading starts as soon as it is ready — usually within a minute."),
        h("div", { class: "reading-steps" },
          h("div", null, h("span", { class: "ok" }, "✓"), `Image prepared · stroke detected: ${source.style === "hand_drawn" ? "hand-drawn" : "clean"}`),
          h("div", null, loaded ? h("span", { class: "ok" }, "✓") : h("span", { class: "spinner small" }), loaded ? `Model loaded · ${model.backend}` : "Loading the model"),
          h("div", null, loaded ? h("span", { class: "spinner small" }) : h("span", { class: "todo" }), "Recognising the structure"),
          h("div", null, h("span", { class: "todo" }), "Deciding whether it deserves a look")),
        h("div", { class: "bar" }, h("div"))));
  }

  // -- giving a structure ------------------------------------------------------------------------------------
  async function upload(url, file) {
    const form = new FormData();
    form.append("file", file);
    return api("POST", url, undefined, form);
  }

  async function uploadPicture(file) {
    state.busy.upload = true;
    state.errors.source = null;
    render();
    try {
      const source = await upload("/api/source/image", file);
      state.sources.picture = source;
      state.result = null;
      state.notThis = false;
      if (!source.read) {
        state.busy.upload = false;
        await readPicture(source);
        return;
      }
      await refreshCheck(false);
    } catch (err) {
      state.errors.source = err.message;
    }
    state.busy.upload = false;
    render();
  }

  async function readPicture(source) {
    state.busy.reading = true;
    state.errors.source = null;
    render();
    try {
      const read = await api("POST", `/api/source/${encodeURIComponent(source.key)}/read`);
      state.sources.picture = read;
      state.busy.reading = false;
      await refreshCheck(false);
    } catch (err) {
      state.errors.source = err.message;
    }
    state.busy.reading = false;
    render();
  }

  async function uploadFile(file) {
    state.busy.upload = true;
    state.errors.source = null;
    render();
    try {
      const answer = await upload("/api/source/file", file);
      state.result = null;
      if (answer.route === "organometallic") {
        state.sources.file = null;
        state.om = { key: answer.key, shown_name: answer.shown_name, state: { substituents: {}, arrangement: 0, oxidation: {}, charge: null }, data: null, multiplicity: "", error: null, drafts: {} };
        state.busy.upload = false;
        await refreshComplex();
        return;
      }
      state.om = null;
      state.sources.file = answer;
      await refreshCheck(false);
    } catch (err) {
      state.errors.source = err.message;
    }
    state.busy.upload = false;
    render();
  }

  async function submitSmiles() {
    const text = state.smilesDraft.trim();
    const current = state.sources.smiles;
    if (!text || (current && current.text === text)) return;
    state.errors.source = null;
    try {
      const source = await api("POST", "/api/source/smiles", { smiles: text });
      source.text = text;
      state.sources.smiles = source;
      state.result = null;
      await refreshCheck(false);
    } catch (err) {
      state.errors.source = err.message;
    }
    render();
  }

  // -- step 2: check (organic) ------------------------------------------------------------------------------------
  function overrideQuery() {
    const a = state.advanced;
    const params = new URLSearchParams();
    if (a.charge !== "") params.set("charge", a.charge);
    if (a.mult !== "") params.set("multiplicity", a.mult);
    if (a.keepAll) params.set("keep_all_fragments", "true");
    return params.toString();
  }

  async function refreshCheck(draw = true) {
    const source = currentSource();
    if (!source || !source.key || (state.kind === "picture" && !source.read)) { if (draw) render(); return; }
    try {
      const query = overrideQuery();
      state.checks[source.key] = await api("GET", `/api/check/${encodeURIComponent(source.key)}${query ? `?${query}` : ""}`);
      state.errors.check = null;
    } catch (err) {
      state.errors.check = err.message;
    }
    if (draw) render();
  }

  function renderOrganic() {
    const source = currentSource();
    const check = currentCheck();
    if (!source || !check) {
      return state.errors.check ? [section("Step 2 · Check the structure", strip("error", state.errors.check))] : [];
    }
    const out = [];
    const gated = check.molecule === null || check.gate.needed;
    out.push(gated ? renderGate(check) : renderSummary(check));
    const passed = check.molecule !== null && !check.errors && (!check.gate.needed || check.gate.confirmed);
    if (check.molecule !== null && check.errors) {
      out.push(section("Step 3 · Output", h("div", { class: "caption" }, "Resolve the errors above before generating an input.")));
    } else if (!passed) {
      out.push(h("div", { class: "caption section" }, "Step 3 · Output unlocks once the structure is confirmed."));
    } else {
      out.push(renderOutput({ formats: ["gaussian", "orca", "xyz", "sdf"], suggested: check.molecule.name, onGenerate: generateOrganic }));
    }
    return out;
  }

  function section(eyebrow, ...children) {
    return h("section", { class: "section stack" }, h("div", { class: "eyebrow" }, eyebrow), ...children);
  }

  function strip(kind, content) {
    return h("div", { class: `strip ${kind}`, role: kind === "error" ? "alert" : null }, h("div", null, content));
  }

  function issueList(issues, skipCodes = []) {
    const shown = (issues || []).filter((i) => !skipCodes.includes(i.code));
    const warnings = shown.filter((i) => i.level !== "info");
    const notes = shown.filter((i) => i.level === "info");
    return [
      ...warnings.map((i) => strip(i.level === "error" ? "error" : "warn", i.message)),
      notes.length ? h("details", { class: "issues" }, h("summary", null, `Details (${notes.length})`), notes.map((i) => h("div", null, i.message))) : null,
    ];
  }

  function stereoSentence(molecule) {
    return molecule.stereo.describe;
  }

  function renderSummary(check) {
    const m = check.molecule;
    const why = check.gate.summary;
    const stereoCount = `${m.stereo.centers} stereocentre${m.stereo.centers === 1 ? "" : "s"}`;
    const head = h("div", { class: "row pad" },
      h("span", { class: "pill good" }, "no check needed"),
      h("span", { class: "mono", style: { fontWeight: 600 } }, m.formula),
      state.expanded
        ? h("span", { class: "mono muted grow", style: { overflowWrap: "anywhere" } }, m.smiles)
        : h("span", { class: "grow muted" }, `${why} · charge ${m.charge} · multiplicity ${m.multiplicity} · ${stereoCount}`),
      h("button", { class: "btn link", "aria-expanded": String(state.expanded), onclick: () => { state.expanded = !state.expanded; render(); } },
        state.expanded ? "Hide depiction ▴" : "See the depiction ▾"));
    const panel = h("div", { class: "panel" }, head);
    if (state.expanded) {
      panel.append(h("div", { class: "summary-body" },
        h("div", { class: "pane", style: { borderRight: "1px solid var(--divider)" } },
          h("div", { class: "pane-box" }, h("img", { src: check.depiction_url, alt: `2D depiction of ${m.formula}, atoms numbered as in the input` }))),
        h("div", { class: "summary-facts" },
          h("div", { class: "facts" },
            fact("Charge", m.charge), fact("Multiplicity", m.multiplicity),
            fact("Stereocentres", m.stereo.centers), fact("Double bonds", m.stereo.bonds)),
          h("div", { class: "prose" }, h("strong", { style: { fontWeight: 600 } }, "Stereochemistry: "), stereoSentence(m)),
          h("div", { class: "small mono" }, `InChIKey ${m.inchikey || "(unavailable)"}`),
          correctionField(check, false))));
    }
    return section("Step 2 · Structure", panel, ...issueList(check.issues));
  }

  const fact = (label, value) => h("div", null, h("div", { class: "fact-label" }, label), h("div", { class: "fact-value" }, String(value)));

  function renderGate(check) {
    const m = check.molecule;
    const panel = h("div", { class: "panel" });
    if (m && check.gate.needed && !check.gate.confirmed) {
      panel.append(h("div", { class: "strip warn" },
        h("span", { class: "pill warn" }, "check"),
        h("div", null, check.gate.reasons.length === 1 ? check.gate.reasons[0]
          : h("ul", null, check.gate.reasons.map((r) => h("li", null, r))))));
    } else if (m && check.gate.confirmed) {
      panel.append(h("div", { class: "strip good", style: { borderRadius: 0, borderWidth: "0 0 1px" } },
        h("span", { class: "pill good" }, "confirmed"), h("div", null, "You confirmed this structure. Correcting the SMILES asks again.")));
    }
    panel.append(h("div", { class: "panes" },
      h("div", { class: "pane" }, h("div", { class: "pane-title" }, "What you gave it"),
        h("div", { class: "pane-box" }, check.preview_url ? h("img", { src: check.preview_url, alt: "The photo you gave" }) : h("span", { class: "muted small" }, "No picture for this structure."))),
      h("div", { class: "pane" }, h("div", { class: "pane-title" }, "What m2i understood"),
        h("div", { class: "pane-box" }, m ? h("img", { src: check.depiction_url, alt: `2D depiction of ${m.formula}, atoms numbered as in the input` })
          : h("div", { class: "strip error", style: { margin: "12px" } }, `This structure cannot be used: ${check.error}`)))));
    if (m) {
      const centres = m.stereo.describe;
      panel.append(h("div", { class: "metrics" },
        metric("Formula", m.formula),
        metric("Charge / mult", `${m.charge} / ${m.multiplicity}`),
        metric("Stereocentres", m.stereo.centers ? `${m.stereo.centers}` : "0", centres),
        metric("Source", `${check.source.backend}${check.source.confidence !== null && check.source.confidence !== undefined ? ` · confidence ${check.source.confidence.toFixed(2)}` : ""}`, null, true)));
    }
    panel.append(correctionField(check, true));
    const said = new Set(check.gate.reasons);
    const rest = (check.issues || []).filter((i) => !said.has(i.message));
    return section("Step 2 · Check the structure", panel, ...issueList(rest, ["recognition.low_confidence"]));
  }

  function metric(label, value, title, wide) {
    return h("div", { class: `metric ${wide ? "wide" : ""}`, title: title || null },
      h("div", { class: "metric-label" }, label), h("div", { class: `metric-value ${wide ? "small" : ""}` }, value));
  }

  function correctionField(check, withButtons) {
    const key = check.key;
    const shown = check.source.corrected ? check.source.correction : check.source.display_smiles;
    const draftKey = `correction:${key}`;
    if (state[draftKey] === undefined) state[draftKey] = shown;
    const input = h("input", {
      class: "field", value: state[draftKey], spellcheck: false, autocomplete: "off", "data-focus": draftKey,
      "aria-label": "SMILES - edit it to correct the reading",
      oninput: (e) => { state[draftKey] = e.target.value; },
      onkeydown: (e) => { if (e.key === "Enter") applyCorrection(check, e.target.value); },
      onchange: (e) => applyCorrection(check, e.target.value),
    });
    const m = check.molecule;
    const buttons = withButtons && m && check.gate.needed && !check.gate.confirmed && !check.errors
      ? h("div", { class: "buttons" },
        h("button", { class: "btn", onclick: () => { state.notThis = true; render(); const el = document.querySelector(`[data-focus="${CSS.escape(draftKey)}"]`); if (el) { el.focus(); el.select(); } } }, "Not this molecule"),
        h("button", { class: "btn primary", onclick: () => confirmStructure(check, true) }, "This is the molecule I drew"))
      : null;
    return h("div", { class: "correction" },
      h("div", { class: "grow" },
        h("label", { class: "label" }, withButtons ? "Wrong? Correct the SMILES and everything downstream follows" : "SMILES — edit it to correct the structure"),
        input,
        state.notThis && withButtons ? h("div", { class: "small", style: { marginTop: "6px" } }, "Type the right SMILES above and press Enter, or give the structure another way: the ChemDraw file, or a clearer photo.") : null,
        check.source.corrected ? h("button", { class: "btn link", style: { marginTop: "6px" }, onclick: () => undoCorrection(check) }, "Back to the original reading") : null),
      buttons);
  }

  async function applyCorrection(check, text) {
    text = text.trim();
    const current = check.source.corrected ? check.source.correction : check.source.display_smiles;
    if (!text || text === current) return;
    try {
      await api("POST", `/api/check/${encodeURIComponent(check.key)}/correct`, { smiles: text });
      state.result = null;
      state.notThis = false;
      await refreshCheck(false);
    } catch (err) {
      state.errors.check = err.message;
    }
    render();
  }

  async function undoCorrection(check) {
    await api("DELETE", `/api/check/${encodeURIComponent(check.key)}/correct`).catch(() => {});
    delete state[`correction:${check.key}`];
    state.result = null;
    await refreshCheck();
  }

  async function confirmStructure(check, confirmed) {
    try {
      await api("POST", `/api/check/${encodeURIComponent(check.key)}/confirm`, { smiles: check.molecule.smiles, confirmed });
      await refreshCheck(false);
    } catch (err) {
      state.errors.check = err.message;
    }
    render();
  }

  // -- step 3: output ---------------------------------------------------------------------------------------------
  async function loadRecipes(program) {
    if (state.output.recipes[program]) return;
    try {
      const data = await api("GET", `/api/profiles?program=${program}`);
      state.output.recipes[program] = data;
      if (data.recipes.length && !state.output.recipe[program]) selectRecipe(program, data.recipes[0].name);
    } catch (err) {
      state.output.problem = err.message;
    }
  }

  function selectRecipe(program, name) {
    const recipe = state.output.recipes[program].recipes.find((r) => r.name === name);
    state.output.recipe[program] = name;
    state.output.fields[program] = {
      method: recipe.method, basis: recipe.basis, dispersion: recipe.dispersion,
      solvent: recipe.solvent || "", nproc: recipe.nproc, mem: recipe.mem,
    };
  }

  function settings() {
    const fmt = state.output.format;
    const fields = state.output.fields[fmt] || {};
    return { format: fmt, recipe: state.output.recipe[fmt] || null, ...fields, nproc: fields.nproc ? Number(fields.nproc) : null };
  }

  let validateTimer = null;
  function scheduleValidate() {
    clearTimeout(validateTimer);
    validateTimer = setTimeout(async () => {
      try {
        const answer = await api("POST", "/api/profile/validate", settings());
        state.output.problem = null;
        state.output.jobs = answer.jobs;
      } catch (err) {
        state.output.problem = err.message;
      }
      const foot = $("output-foot");
      if (foot) foot.replaceWith(outputFoot(state.outputContext));
    }, 350);
  }

  function renderOutput(context) {
    state.outputContext = context;
    const fmt = context.formats.includes(state.output.format) ? state.output.format : context.formats[0];
    state.output.format = fmt;
    const program = state.output.recipes[fmt];
    const panel = h("div", { class: "panel panel-pad" },
      h("div", { class: "formats", role: "radiogroup", "aria-label": "Format" },
        context.formats.map((f) => h("button", {
          class: `format ${f === fmt ? "on" : ""}`, role: "radio", "aria-checked": String(f === fmt),
          onclick: async () => { state.output.format = f; await loadRecipes(f); state.output.problem = null; render(); scheduleValidate(); },
        }, `${FORMATS[f][1]} `, h("span", null, FORMATS[f][0])))));

    if (["gaussian", "orca"].includes(fmt) && program) {
      const fields = state.output.fields[fmt] || {};
      const recipes = program.recipes;
      const current = recipes.find((r) => r.name === state.output.recipe[fmt]);
      const set = (key) => (e) => { fields[key] = e.target.value; scheduleValidate(); };
      panel.append(
        h("div", { style: { marginBottom: "16px" } },
          h("label", { class: "label", for: "recipe" }, "Recipe"),
          h("select", { class: "field tall", id: "recipe", style: { maxWidth: "520px" }, onchange: (e) => { selectRecipe(fmt, e.target.value); render(); scheduleValidate(); } },
            recipes.map((r) => h("option", { value: r.name, selected: current && r.name === current.name }, r.name))),
          current && current.description ? h("div", { class: "small", style: { marginTop: "6px" } }, current.description) : null),
        h("div", { class: "grid3" },
          textField("Method", fields.method, set("method"), `method-${fmt}`),
          textField("Basis set", fields.basis, set("basis"), `basis-${fmt}`),
          h("div", null, h("label", { class: "label" }, "Dispersion"),
            h("select", { class: "field", onchange: set("dispersion") },
              program.dispersion.map((d) => h("option", { value: d, selected: d === (fields.dispersion || "none") }, d)))),
          textField("Solvent (blank = gas phase)", fields.solvent, set("solvent"), `solvent-${fmt}`, "gas phase"),
          textField("Cores", fields.nproc, set("nproc"), `nproc-${fmt}`, null, "number"),
          textField("Memory", fields.mem, set("mem"), `mem-${fmt}`)));
    } else if (fmt === "xyz") {
      panel.append(h("div", { class: "caption" }, "Geometry only: element symbols and Cartesian coordinates, readable by practically every program."));
    } else if (fmt === "sdf") {
      panel.append(h("div", { class: "caption" }, "Geometry only, but it keeps bonds, charges and stereochemistry: the format to hand to another cheminformatics tool."));
    }
    panel.append(nameField(context));
    panel.append(outputFoot(context));
    const out = [section(context.eyebrow || "Step 3 · Output", panel)];
    if (context.blocked) out.push(strip("warn", context.blocked));
    if (state.errors.output) out.push(strip("error", state.errors.output));
    if (state.result && resultCurrent()) out.push(renderResult(state.result));
    return h("div", { class: "section" }, ...out);
  }

  function textField(label, value, oninput, focus, placeholder, type) {
    return h("div", null, h("label", { class: "label" }, label),
      h("input", { class: "field", type: type || "text", min: type === "number" ? 1 : null, value: value === undefined || value === null ? "" : value, placeholder, "data-focus": focus, oninput, spellcheck: false }));
  }

  // The written files take this name; blank leaves it to m2i (the formula, the
  // molecule's own name, the CIF it came from). Cleaned here as the server does.
  function cleanName(text) {
    return (text || "").trim().replace(/[^A-Za-z0-9_-]+/g, "_").replace(/^_+|_+$/g, "").slice(0, 60);
  }

  function nameField(context) {
    const suggested = context.suggested || "molecule";
    const note = h("div", { class: "small", style: { marginTop: "6px" } }, nameHint(suggested));
    return h("div", { style: { marginTop: "16px", maxWidth: "520px" } },
      textField("File name", state.output.name, (e) => {
        state.output.name = e.target.value;
        note.textContent = nameHint(suggested);
      }, "output-name", suggested),
      note);
  }

  function nameHint(suggested) {
    const stem = cleanName(state.output.name) || suggested;
    const several = showAdvanced() && Number(state.advanced.keep) > 1;
    const files = `${stem}${several ? "_c01" : ""}${FORMATS[state.output.format][0]}`;
    return cleanName(state.output.name) ? `Written as ${files}` : `Written as ${files}, unless you name it yourself`;
  }

  function outputFoot(context) {
    const fmt = state.output.format;
    const busy = state.busy.generate;
    return h("div", { class: "output-foot", id: "output-foot" },
      state.output.problem ? h("div", { class: "strip error grow", style: { padding: "8px 12px" } }, state.output.problem)
        : h("div", { class: "mono muted", style: { fontSize: "12px" } },
          ["gaussian", "orca"].includes(fmt) ? `jobs from the recipe: ${state.output.jobs || jobsOf(fmt)}` : "geometry only"),
      h("button", { class: "btn primary big", disabled: Boolean(busy || state.output.problem || context.blocked), onclick: context.onGenerate },
        busy ? [h("span", { class: "spinner small", style: { display: "inline-block", marginRight: "8px", borderColor: "#c7d2fe", borderTopColor: "#fff" } }), "Generating…"] : `Generate ${FORMAT_LABEL[fmt]}`));
  }

  function jobsOf(fmt) {
    const program = state.output.recipes[fmt];
    const recipe = program && program.recipes.find((r) => r.name === state.output.recipe[fmt]);
    return recipe ? recipe.jobs.join(" ") || "single point" : "";
  }

  function signature() {
    const source = currentSource();
    const check = currentCheck();
    const base = { kind: state.kind, settings: settings(), name: cleanName(state.output.name) };
    if (state.kind === "file" && state.om) return JSON.stringify({ ...base, om: state.om.key, st: state.om.state, mult: state.om.multiplicity });
    if (state.kind === "cif" && state.cif) return JSON.stringify({ ...base, cif: state.cif.summary && state.cif.summary.key, st: state.cif.state, mult: state.cif.multiplicity });
    return JSON.stringify({ ...base, key: source && source.key, smiles: check && check.molecule && check.molecule.smiles, advanced: state.advanced });
  }

  function resultCurrent() {
    return state.result && state.result.signature === signature();
  }

  async function runGenerate(url, body) {
    state.busy.generate = true;
    state.errors.output = null;
    render();
    const sig = signature();
    try {
      const result = await api("POST", url, body);
      result.signature = sig;
      state.result = result;
    } catch (err) {
      state.errors.output = err.message;
    }
    state.busy.generate = false;
    render();
    const el = $("result");
    if (el) el.scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function generateOrganic() {
    const source = currentSource();
    const a = state.advanced;
    return runGenerate("/api/generate", {
      key: source.key,
      settings: settings(),
      name: cleanName(state.output.name) || null,
      conformers: { n_confs: a.nconfs === "" ? null : Number(a.nconfs), keep: Number(a.keep) || 1, seed: Number(a.seed) || 61453, force_field: a.forceField },
      overrides: { charge: a.charge === "" ? null : Number(a.charge), multiplicity: a.mult === "" ? null : Number(a.mult), keep_all_fragments: a.keepAll },
    });
  }

  function renderResult(result) {
    const files = result.files;
    const first = files[0];
    const conformer = first && first.conformer;
    const sub = conformer
      ? `conformer 1 of ${files.length} · ${conformer.force_field} · geometry reused if you switch format`
      : "experimental or built geometry, as shown above";
    const wrap = h("div", { class: "section stack", id: "result" },
      h("div", { class: "result-bar" },
        h("span", { class: "tick" }, "✓"),
        h("div", { class: "grow" }, h("div", { class: "result-title" }, `${files.length} ${result.format_label} file${files.length === 1 ? "" : "s"} ready`), h("div", { class: "result-sub" }, sub)),
        first ? h("a", { class: "btn success big", href: first.download_url, download: first.name }, `Download ${first.name}`) : null),
      ...result.issues.map((i) => strip(i.level === "error" ? "error" : "warn", i.message)));
    files.forEach((file, index) => {
      const open = state.previewOpen === undefined ? index === 0 : state.previewOpen === file.name;
      const c = file.conformer;
      wrap.append(h("div", { class: "panel" },
        h("div", { class: "preview-head" },
          h("span", { class: "mono" }, file.name),
          c ? h("span", { class: "small mono" }, `conformer ${c.index + 1}${c.relative_energy !== null && c.relative_energy !== undefined ? `, +${c.relative_energy.toFixed(2)} kcal/mol (${c.force_field})` : ""}`) : null,
          h("span", { class: "grow" }),
          index > 0 ? h("a", { class: "btn ghost", href: file.download_url, download: file.name }, "Download") : null,
          h("button", { class: "btn link", "aria-expanded": String(open), onclick: () => { state.previewOpen = open ? null : file.name; render(); } }, open ? "Collapse ▴" : "Preview ▾")),
        open ? h("pre", { class: "preview" }, file.preview) : null));
    });
    if (result.extras.length) {
      wrap.append(h("div", { class: "extras" }, result.extras.map((extra) => h("div", { class: "extra" },
        h("div", { class: "file-badge" }, extra.kind === "provenance" ? "JSON" : "PNG"),
        h("div", { class: "grow" }, h("div", { class: "file-name" }, extra.name),
          h("div", { class: "file-meta" }, extra.kind === "provenance" ? "versions · source · confidence · profile · every warning and decision" : "atoms numbered as in the input")),
        h("a", { class: "btn ghost", href: extra.download_url, download: extra.name }, "Download")))));
    }
    return wrap;
  }

  // -- 3D view -------------------------------------------------------------------------------------------------------
  const viewers = new Map();
  function viewer(payload) {
    if (viewers.has(payload.xyz)) return viewers.get(payload.xyz);
    const el = h("div", { class: "viewer", role: "img", "aria-label": "3D view of the structure; drag to rotate, scroll to zoom" });
    if (typeof window.$3Dmol === "undefined") {
      el.append(h("div", { class: "viewer-note" }, "3D view unavailable"));
      return el;
    }
    viewers.clear(); // one structure at a time
    viewers.set(payload.xyz, el);
    requestAnimationFrame(() => {
      const view = window.$3Dmol.createViewer(el, { backgroundColor: "#14161b" });
      view.addModel(payload.xyz, "xyz");
      view.setStyle({}, { stick: { radius: 0.14 }, sphere: { scale: 0.22 } });
      (payload.labels || []).forEach((label) => view.addLabel(label.text, {
        position: { x: label.x, y: label.y, z: label.z }, fontSize: 12, fontColor: "white", backgroundColor: "#4338ca", backgroundOpacity: 0.8,
      }));
      view.zoomTo();
      view.render();
      new ResizeObserver(() => { view.resize(); view.render(); }).observe(el);
    });
    return el;
  }

  function coordinationCard(centres, extra) {
    return h("div", { class: "panel coord" },
      centres.length ? centres.map((centre) => h("div", { style: { marginBottom: "10px" } },
        h("div", { class: "coord-title" }, `${centre.label} — ${centre.geometry}`, centre.detail ? h("span", { class: "muted", style: { fontWeight: 400 } }, ` (${centre.detail})`) : null),
        h("div", { class: "table-scroll" }, h("table", { class: "table" },
          h("thead", null, h("tr", null, h("th", null, "Donor"), h("th", null, "Distance (Å)"), h("th", null, "Type"))),
          h("tbody", null, centre.donors.map((d) => h("tr", null, h("td", { class: "mono" }, d.label), h("td", { class: "mono" }, d.distance.toFixed(3)), h("td", null, d.kind)))))),
        centre.trans_pairs.map((t) => h("div", { class: "small" }, `trans: ${t.a} / ${t.b} (${t.angle}°)`)),
        centre.isomers.map((iso) => h("div", { class: "small", style: { fontWeight: 600, color: "var(--body)" } }, iso)),
        centre.helicity ? h("div", { class: "small" }, `${centre.helicity} at ${centre.label} (${centre.helicity_detail})`) : null))
        : h("div", { class: "small" }, "No metal in this species."),
      extra);
  }

  // -- electronic state (complexes and crystals) --------------------------------------------------------------------
  function electronicSection(eyebrow, electronic, holder, refresh, intro) {
    const panel = h("div", { class: "panel panel-pad" });
    if (intro) panel.append(h("div", { class: "prose", style: { marginBottom: "12px" } }, intro));
    const grid = h("div", { class: "grid3" });
    electronic.metals.forEach((metal) => {
      grid.append(h("div", null, h("label", { class: "label" }, `Oxidation state of ${metal}`),
        h("input", { class: "field", type: "number", step: 1, placeholder: "e.g. 2", value: holder.state.oxidation[metal] ?? "", "data-focus": `ox-${metal}`,
          onchange: (e) => { if (e.target.value === "") delete holder.state.oxidation[metal]; else holder.state.oxidation[metal] = Number(e.target.value); refresh(); } })));
    });
    grid.append(h("div", null, h("label", { class: "label" }, "Charge"),
      h("input", { class: "field", type: "number", step: 1, value: holder.state.charge ?? electronic.charge ?? "", placeholder: "required", "data-focus": "charge",
        onchange: (e) => { holder.state.charge = e.target.value === "" ? null : Number(e.target.value); refresh(); } })));
    if (electronic.options.length) {
      grid.append(h("div", null, h("label", { class: "label" }, "Spin state"),
        h("select", { class: "field", onchange: (e) => { holder.multiplicity = e.target.value; render(); } },
          h("option", { value: "", selected: holder.multiplicity === "" }, "choose"),
          electronic.options.map((o) => h("option", { value: o.multiplicity, selected: String(o.multiplicity) === String(holder.multiplicity) }, `${o.multiplicity} — ${o.label} (${o.unpaired} unpaired)`)))));
    } else {
      grid.append(h("div", null, h("label", { class: "label" }, "Multiplicity"),
        h("input", { class: "field", type: "number", min: 1, step: 1, value: holder.multiplicity, placeholder: electronic.lowest ? `lowest allowed: ${electronic.lowest}` : "", "data-focus": "mult",
          oninput: (e) => { holder.multiplicity = e.target.value; },
          onchange: () => render() })));
    }
    panel.append(grid,
      h("div", { class: "small", style: { marginTop: "10px" } }, `Suggestion: ${electronic.charge_reason}.`),
      ...electronic.hints.map((line) => h("div", { class: "small" }, line)));
    return section(eyebrow, panel);
  }

  // -- metal complexes (2e) --------------------------------------------------------------------------------------------
  let complexTimer = null;
  async function refreshComplex() {
    const om = state.om;
    om.busy = true;
    render();
    try {
      om.data = await api("POST", `/api/om/${encodeURIComponent(om.key)}/state`, om.state);
      om.error = null;
      om.state.arrangement = om.data.arrangement;
    } catch (err) {
      om.error = err.message;
    }
    om.busy = false;
    render();
  }

  function complexSteps() {
    const om = state.om;
    const d = om.data;
    const result = state.result && resultCurrent();
    const chosen = d && d.arrangements[d.arrangement];
    return [
      { name: "Structure", sub: d ? `${d.drawing.metal} · ${d.drawing.donors.length} donors` : "reading the drawing…", status: d ? "done" : "active" },
      { name: "Check", sub: d ? (d.ambiguous ? "choose the isomer" : chosen ? chosen.label.split(";")[0] : "") : null, subClass: d && d.ambiguous ? "warn" : "", status: d ? "done" : "pending" },
      { name: "Electronic state", sub: om.multiplicity ? `mult ${om.multiplicity}` : d ? "charge and spin" : null, status: om.multiplicity ? "done" : d ? "active" : "pending" },
      { name: "Output", sub: result ? `${state.result.format} · ${state.result.jobs}` : null, status: result ? "done" : om.multiplicity ? "active" : "pending" },
    ];
  }

  function renderComplex() {
    const om = state.om;
    const d = om.data;
    if (!d) return [om.error ? strip("error", om.error) : h("div", { class: "row section" }, h("span", { class: "spinner" }), "Reading the drawing and building the complex…")];
    const out = [];
    const summary = h("div", { class: "panel", style: { padding: "13px 16px" } },
      h("div", { style: { font: "600 13px var(--sans)" } }, `Metal complex: ${d.drawing.metal} with ${d.drawing.donors.length} donor atoms — ${d.drawing.donors.map((n) => n.split(" of ")[0]).join(", ")}`),
      h("div", { class: "small mono" }, `read from ${d.shown_name} · bonds to the metal taken as drawn`));
    out.push(h("div", { class: "section stack" }, summary, ...d.drawing.notes.map((n) => n.level === "warning" ? strip("warn", n.message) : h("div", { class: "small" }, n.message))));

    if (d.drawing.gaps.length) {
      const rows = d.drawing.gaps.map((gap) => {
        const draft = om.drafts[gap.atom] ?? om.state.substituents[gap.atom] ?? "";
        return h("tr", null, h("td", { class: "mono" }, gap.atom), h("td", { class: "mono" }, gap.bonds_missing),
          h("td", null, h("input", { class: "field", value: draft, placeholder: "blank leaves hydrogens", "data-focus": `gap-${gap.atom}`, spellcheck: false,
            oninput: (e) => { om.drafts[gap.atom] = e.target.value; },
            onkeydown: (e) => { if (e.key === "Enter") e.target.blur(); },
            onchange: (e) => { om.state.substituents[gap.atom] = e.target.value.trim(); if (!om.state.substituents[gap.atom]) delete om.state.substituents[gap.atom]; om.state.arrangement = 0; refreshComplex(); } })));
      });
      out.push(section("What the drawing leaves out",
        h("div", { class: "caption", style: { marginBottom: "10px" } }, "Write the real groups — iPr2, Ph2, Cy2, Me, OMe, CH2OH, or SMILES — one per missing bond. Left empty, the valence is filled with hydrogen, which is hardly ever what a scheme means."),
        h("div", { class: "panel table-scroll" }, h("table", { class: "table" },
          h("thead", null, h("tr", null, h("th", null, "Atom"), h("th", null, "Bonds missing"), h("th", { style: { width: "55%" } }, "Groups"))),
          h("tbody", null, rows)))));
    }

    const check = section("Step 2 · Check the structure");
    if (om.error) check.append(strip("error", om.error));
    if (d.ambiguous) {
      check.append(strip("warn", d.mirror_only
        ? "The first two are mirror images, and the drawing does not tell them apart: only wedges and hashes on the bonds to the metal do. Choose the enantiomer."
        : "The drawing fits the first two about equally well: choose one."));
      check.lastChild.style.marginBottom = "10px";
    }
    check.append(h("div", { class: "choices", role: "radiogroup", "aria-label": "Arrangement of the ligands" },
      d.arrangements.map((a) => h("label", { class: `choice ${a.index === d.arrangement ? "selected" : ""}` },
        h("input", { type: "radio", name: "arrangement", checked: a.index === d.arrangement, onchange: () => { om.state.arrangement = a.index; refreshComplex(); } }),
        h("span", null, `${a.index + 1}. ${a.label}${a.mirror_of !== null && a.mirror_of !== undefined ? ` — mirror image of ${a.mirror_of + 1}` : ""}`),
        h("span", { class: "fit" }, `fit ${a.fit.toFixed(2)}`)))));
    check.append(h("div", { class: "split", style: { marginTop: "14px" } },
      h("div", null, viewer(d.built.viewer), h("div", { class: "small", style: { marginTop: "6px" } }, `${d.built.formula}, built by distance geometry and pre-optimised with UFF: a starting geometry for your optimisation.`)),
      coordinationCard(d.built.coordination)));
    d.built.notes.forEach((note) => check.append(strip("warn", note)));
    out.push(check);

    out.push(electronicSection("Step 3 · Electronic state", d.electronic, om, refreshComplex));
    const blocked = !om.multiplicity ? "Choose the spin state to continue." : null;
    out.push(renderOutput({ eyebrow: "Step 4 · Output", formats: ["gaussian", "orca", "xyz"], blocked, suggested: d.built.name, onGenerate: generateComplex }));
    return out;
  }

  function generateComplex() {
    const om = state.om;
    return runGenerate(`/api/om/${encodeURIComponent(om.key)}/generate`, {
      ...om.state, charge: om.state.charge ?? om.data.electronic.charge, multiplicity: Number(om.multiplicity), settings: settings(), name: cleanName(state.output.name) || null,
    });
  }

  // -- crystals (2f) --------------------------------------------------------------------------------------------------------
  async function uploadCif(file, normalise = true) {
    state.busy.upload = true;
    state.errors.source = null;
    render();
    try {
      const summary = await upload(`/api/source/cif?normalise=${normalise}`, file);
      const species = summary.species.findIndex((s) => s.metals.length) >= 0 ? summary.species.findIndex((s) => s.metals.length) : 0;
      state.cif = { file, summary, data: null, error: null, multiplicity: "",
        state: { species, counts: {}, add_hydrogens: null, allow_missing: false, priorities: {}, mirror: false, oxidation: {}, charge: null } };
      state.result = null;
      state.busy.upload = false;
      await refreshCrystal();
      return;
    } catch (err) {
      state.errors.source = err.message;
    }
    state.busy.upload = false;
    render();
  }

  async function refreshCrystal() {
    const cif = state.cif;
    try {
      cif.data = await api("POST", `/api/crystal/${encodeURIComponent(cif.summary.key)}/state`, cif.state);
      cif.error = null;
    } catch (err) {
      cif.error = err.message;
    }
    render();
  }

  function crystalSteps() {
    const cif = state.cif;
    const d = cif.data;
    const result = state.result && resultCurrent();
    const hydrogensPending = d && d.hydrogens.needed && !d.hydrogens.add && !cif.state.allow_missing;
    return [
      { name: "Structure", sub: `cif · ${cif.summary.species.length} species`, status: "done" },
      { name: "Check", sub: hydrogensPending ? "hydrogens missing" : d ? "complete" : null, subClass: hydrogensPending ? "warn" : "good", status: d && !hydrogensPending ? "done" : "active" },
      { name: "Electronic state", sub: cif.multiplicity ? `mult ${cif.multiplicity}` : null, status: cif.multiplicity ? "done" : d && !hydrogensPending ? "active" : "pending" },
      { name: "Output", sub: result ? `${state.result.format} · ${state.result.jobs}` : null, status: result ? "done" : cif.multiplicity ? "active" : "pending" },
    ];
  }

  function renderCrystalFile() {
    const s = state.cif.summary;
    const [a, b, c, alpha, beta, gamma] = s.cell;
    const cell = `${s.source} · ${s.spacegroup}${s.z ? ` · Z = ${s.z}` : ""} · a ${a}, b ${b}, c ${c} Å, α ${alpha}°, β ${beta}°, γ ${gamma}°`;
    return h("div", null,
      h("div", { class: "file-row" },
        h("div", { class: "file-badge" }, "CIF"),
        h("div", { class: "grow" }, h("div", { class: "file-name" }, s.shown_name), h("div", { class: "file-meta" }, cell)),
        h("label", { class: "check" }, h("input", { type: "checkbox", checked: s.normalise, onchange: (e) => uploadCif(state.cif.file, e.target.checked) }), "Extend C–H, N–H and O–H to standard lengths"),
        h("button", { class: "btn ghost", onclick: () => { state.cif = null; state.result = null; render(); } }, "Replace")),
      h("div", { class: "small" }, "X-ray places hydrogens about 0.1 Å too close to their atom. Hydrides and other M–H are never moved."));
  }

  function renderCrystal() {
    const cif = state.cif;
    if (!cif || !cif.summary) return [];
    const s = cif.summary;
    const d = cif.data;
    const out = [];
    out.push(h("div", { class: "section" },
      h("div", { class: "panel table-scroll" }, h("table", { class: "table" },
        h("thead", null, h("tr", null, h("th", null, "Species"), h("th", null, "Copies in the cell"), h("th", null, "Atoms"), h("th", null, "Metal"), h("th", null, "Calculate"))),
        h("tbody", null, s.species.map((sp) => h("tr", { class: sp.index === cif.state.species ? "selected" : "" },
          h("td", { class: "mono" }, sp.formula), h("td", { class: "mono" }, sp.copies), h("td", { class: "mono" }, sp.atoms),
          h("td", { class: "mono" }, sp.metals.join(", ") || "—"),
          h("td", null, h("input", { type: "radio", name: "species", checked: sp.index === cif.state.species, "aria-label": `Calculate ${sp.formula}`,
            onchange: () => { cif.state = { ...cif.state, species: sp.index, counts: {}, add_hydrogens: null, priorities: {}, oxidation: {}, charge: null }; cif.multiplicity = ""; refreshCrystal(); } }))))))),
      ...s.extended.map((formula) => h("div", { class: "small", style: { marginTop: "6px" } }, `${formula}: an extended network with no discrete molecule; it cannot be calculated as one.`)),
      ...issueList(s.issues)));
    if (cif.error) out.push(strip("error", cif.error));
    if (!d) return out;

    const check = section("Step 2 · Check the structure");
    const hyd = d.hydrogens;
    if (hyd.needed) {
      const status = hyd.target === null
        ? strip("warn", `${hyd.total} hydrogens proposed. The formula does not say how many this species needs, so check the list.`)
        : hyd.matches ? strip("good", `${hyd.total} hydrogens, as the formula requires.`)
          : strip("error", `${hyd.total} hydrogens proposed, but the formula needs ${hyd.target}. Adjust the rows marked 'assumed'.`);
      check.append(h("div", { class: "panel", style: { marginBottom: "14px" } },
        h("div", { class: "coord-title", style: { padding: "12px 14px 0" } }, "Missing hydrogens"),
        h("div", { class: "table-scroll" }, h("table", { class: "table" },
          h("thead", null, h("tr", null, h("th", null, "Atom"), h("th", null, "H to add"), h("th", null, "Why"), h("th", null, "Sure"))),
          h("tbody", null, hyd.sites.map((site) => h("tr", null,
            h("td", { class: "mono" }, site.label),
            h("td", { style: { width: "110px" } }, h("input", { class: "field", type: "number", min: 0, max: 4, step: 1, value: site.count, "aria-label": `Hydrogens to add on ${site.label}`,
              onchange: (e) => { cif.state.counts[site.label] = Number(e.target.value); cif.state.add_hydrogens = null; refreshCrystal(); } })),
            h("td", null, site.reason),
            h("td", null, h("span", { class: site.certain ? "sure-yes" : "sure-no" }, site.certain ? "yes" : "assumed"))))))),
        h("div", { style: { padding: "12px 14px" }, class: "stack" }, status,
          h("label", { class: "check" }, h("input", { type: "checkbox", checked: hyd.add, onchange: (e) => { cif.state.add_hydrogens = e.target.checked; refreshCrystal(); } }), "Add these hydrogens"),
          !hyd.add ? h("label", { class: "check" }, h("input", { type: "checkbox", checked: cif.state.allow_missing, onchange: (e) => { cif.state.allow_missing = e.target.checked; refreshCrystal(); } }), "Write the input without them — I know the hydrogens are missing") : null)));
    } else if (hyd.elsewhere) {
      check.append(h("div", { class: "small", style: { marginBottom: "10px" } }, `The hydrogens missing from the crystal belong to other species; ${d.species.formula} is complete.`));
    }
    const chir = d.chirality;
    const extra = h("div", { class: "stack", style: { marginTop: "8px" } },
      ...chir.questions.map((q) => h("div", null,
        h("div", { class: "small" }, `Ring ${q.ring}: ${q.reason}.`),
        h("div", { class: "row", style: { gap: "16px" } }, h("span", { class: "small" }, `Higher CIP priority in ring ${q.ring}:`),
          q.options.map((option) => h("label", { class: "check" }, h("input", { type: "radio", name: `ring-${q.ring}`, checked: q.chosen === option,
            onchange: () => { cif.state.priorities[q.ring] = option; refreshCrystal(); } }), option))))),
      ...chir.notes.map((note) => h("div", { class: "small" }, note)),
      chir.chiral ? h("label", { class: "check" }, h("input", { type: "checkbox", checked: cif.state.mirror, onchange: (e) => { cif.state.mirror = e.target.checked; render(); } }), "Write the mirror image (the other enantiomer)") : null);
    check.append(h("div", { class: "split" }, viewer(d.viewer), coordinationCard(d.coordination, extra)));
    out.push(check);

    const e = d.electronic;
    const intro = e.needs_answers
      ? `${d.species.formula} contains ${e.metals.join(", ")}. A crystal does not record oxidation states, charge or spin, so these are yours to set.` : null;
    out.push(electronicSection("Step 3 · Electronic state", { ...e, charge: e.charge }, cif, refreshCrystal, intro));
    const blocked = d.blockers[0] || (cif.state.charge === null && e.charge === null ? "Set the charge to continue." : null) || (!cif.multiplicity ? "Choose the spin state to continue." : null);
    out.push(renderOutput({ eyebrow: "Step 4 · Output", formats: ["gaussian", "orca", "xyz"], blocked, suggested: d.species.name, onGenerate: generateCrystal }));
    return out;
  }

  function generateCrystal() {
    const cif = state.cif;
    return runGenerate(`/api/crystal/${encodeURIComponent(cif.summary.key)}/generate`, {
      ...cif.state, charge: cif.state.charge ?? cif.data.electronic.charge, multiplicity: Number(cif.multiplicity), settings: settings(), name: cleanName(state.output.name) || null,
    });
  }

  // -- splash (3a) --------------------------------------------------------------------------------------------------------------
  function renderSplash() {
    const splash = $("splash");
    const health = state.health;
    const model = health && health.model;
    const show = model && (model.state === "loading" || model.state === "waiting") && !state.splashDismissed && nothingGiven();
    splash.hidden = !show;
    document.body.style.overflow = show ? "hidden" : "";
    if (!show) return;
    const enter = (kind) => () => { state.splashDismissed = true; chooseKind(kind); };
    const fill = Math.min(95, Math.round((model.elapsed / model.typical) * 100));
    splash.replaceChildren(h("div", { class: "splash-inner", role: "status", "aria-live": "polite" },
      h("img", { class: "splash-logo", src: "logo.svg", alt: "m2i" }),
      h("h1", null, "Starting m2i"),
      h("p", null, "The server has just started. It loads the recognition model and is ready — usually in under a minute. No need to reload: the page comes in on its own."),
      h("div", { class: "splash-steps" },
        splashStep("Server running", `${health.uptime} s`, "done", true),
        splashStep(`Loading ${model.backend ? model.backend.toUpperCase() : "the recognition model"}`,
          [h("div", { class: "splash-detail" }, `2.4 GB in memory · loaded once, for everyone · ${model.elapsed} s of about ${model.typical} s`),
            h("div", { class: "bar fill" }, h("div", { style: { width: `${fill}%` } }))], "running", true),
        splashStep("Ready to read", null, "pending", false)),
      h("div", { class: "splash-more" },
        h("div", { class: "eyebrow" }, "In the meantime"),
        h("div", { class: "text" }, "Starting from a ", h("strong", null, "SMILES"), ", a ", h("strong", null, "ChemDraw"), " file or a ", h("strong", null, ".cif"), "? None of those need the model — those routes are open already."),
        h("div", { class: "splash-buttons" },
          h("button", { class: "btn dark", onclick: enter("smiles") }, "Start with SMILES"),
          h("button", { class: "btn dark", onclick: enter("file") }, "Upload ChemDraw"),
          h("button", { class: "btn dark", onclick: enter("cif") }, "Upload .cif")))));
  }

  function splashStep(name, detail, status, line) {
    return h("div", { class: `splash-step ${status === "pending" ? "off" : ""}` },
      h("div", { class: "step-track" },
        h("span", { class: `step-dot ${status}` }, status === "done" ? "✓" : ""),
        line ? h("span", { class: "step-line" }) : null),
      h("div", { class: "grow" }, h("div", { class: "splash-name" }, name),
        typeof detail === "string" ? h("div", { class: "splash-detail" }, detail) : detail));
  }

  $("home").addEventListener("click", (event) => {
    event.preventDefault();
    window.location.reload();
  });

  boot();
})();
