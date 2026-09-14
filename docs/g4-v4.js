(() => {
  "use strict";
  const manifest = window.G4V4_MANIFEST;
  const cases = window.G4V4_CASES || {};
  const events = window.G4V4_EVENTS || {};
  if (!manifest) {
    document.getElementById("interpretation").textContent = "Generated data are missing. Run the G4 v4 builder first.";
    return;
  }

  const $ = id => document.getElementById(id);
  const canvas = $("network");
  const ctx = canvas.getContext("2d");
  const plot = $("plot");
  const pctx = plot.getContext("2d");
  const eventCanvas = $("eventCanvas");
  const ectx = eventCanvas.getContext("2d");
  const initial = manifest.geometry.initial;
  const edges = manifest.geometry.edges;
  const fibers = manifest.geometry.fibers;
  const fixed = new Set(manifest.geometry.fixed);
  const domain = manifest.domain;
  const decoded = new Map();
  const eventDecoded = new Map();
  let caseId = "normal";
  let frameIndex = 0;
  let playing = false;
  let timer = null;

  function loadScript(path) {
    return new Promise((resolve, reject) => {
      const script = document.createElement("script");
      script.src = path;
      script.onload = resolve;
      script.onerror = () => reject(new Error(`Could not load ${path}`));
      document.head.appendChild(script);
    });
  }

  async function loadCase(id) {
    if (!cases[id]) await loadScript(manifest.caseFiles[id]);
    if (id === "clutch" && !events.clutch) await loadScript(manifest.eventFiles.clutch);
    return cases[id];
  }

  function unpack(pack, base) {
    const binary = atob(pack.base64);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
    const values = pack.dtype === "i1" ? new Int8Array(bytes.buffer) : new Int16Array(bytes.buffer);
    const [nFrames, nPoints] = pack.shape;
    const frames = new Array(nFrames);
    let cursor = 0;
    for (let f = 0; f < nFrames; f++) {
      const points = new Array(nPoints);
      for (let i = 0; i < nPoints; i++) {
        points[i] = [
          base[i][0] + pack.scale * values[cursor++],
          base[i][1] + pack.scale * values[cursor++],
        ];
      }
      frames[f] = points;
    }
    return frames;
  }

  function positionsFor(id) {
    if (!decoded.has(id)) decoded.set(id, unpack(cases[id].positions, initial));
    return decoded.get(id);
  }

  function eventPositions(event) {
    const key = caseId;
    if (!eventDecoded.has(key)) {
      const base = event.localBeads.map(i => initial[i]);
      eventDecoded.set(key, unpack(event.positions, base));
    }
    return eventDecoded.get(key);
  }

  function xy(point, target = canvas) {
    const pad = 24;
    const scale = (Math.min(target.width, target.height) - 2 * pad) / domain;
    return [target.width / 2 + point[0] * scale, target.height / 2 - point[1] * scale];
  }

  function xyLocal(point, target, center, span) {
    const pad = 28;
    const scale = (Math.min(target.width, target.height) - 2 * pad) / span;
    return [
      target.width / 2 + (point[0] - center[0]) * scale,
      target.height / 2 - (point[1] - center[1]) * scale,
    ];
  }

  function clear(context, target) {
    context.clearRect(0, 0, target.width, target.height);
    context.fillStyle = "#f8f5ed";
    context.fillRect(0, 0, target.width, target.height);
  }

  function line(context, a, b, color, width = 1, alpha = 1) {
    context.save();
    context.globalAlpha = alpha;
    context.strokeStyle = color;
    context.lineWidth = width;
    context.beginPath();
    context.moveTo(a[0], a[1]);
    context.lineTo(b[0], b[1]);
    context.stroke();
    context.restore();
  }

  function arrow(context, origin, vector, color, width = 1.5, multiplier = 1) {
    const magnitude = Math.hypot(vector[0], vector[1]);
    if (!Number.isFinite(magnitude) || magnitude < 1e-12) return;
    const display = Math.min(25, magnitude * multiplier);
    const dx = vector[0] / magnitude * display;
    const dy = -vector[1] / magnitude * display;
    const end = [origin[0] + dx, origin[1] + dy];
    line(context, origin, end, color, width, .9);
    const angle = Math.atan2(dy, dx);
    const size = 4;
    context.fillStyle = color;
    context.beginPath();
    context.moveTo(end[0], end[1]);
    context.lineTo(end[0] - size * Math.cos(angle - .55), end[1] - size * Math.sin(angle - .55));
    context.lineTo(end[0] - size * Math.cos(angle + .55), end[1] - size * Math.sin(angle + .55));
    context.closePath();
    context.fill();
  }

  function materialPoint(positions, link, side) {
    const edgeId = side === 0 ? link[0] : link[2];
    const alpha = side === 0 ? link[1] : link[3];
    const edge = edges[edgeId];
    const a = positions[edge[0]], b = positions[edge[1]];
    return [(1 - alpha) * a[0] + alpha * b[0], (1 - alpha) * a[1] + alpha * b[1]];
  }

  function fiberColor(distance) {
    if (distance === 0) return "#d45e38";
    if (distance === 1) return "#426f97";
    if (distance >= 2) return "#8055a3";
    return "#71808a";
  }

  function drawNetwork() {
    const currentCase = cases[caseId];
    const positions = positionsFor(caseId)[frameIndex];
    const metric = currentCase.metrics[frameIndex] || {};
    const center = currentCase.centers[frameIndex] || [0, 0];
    const distances = currentCase.fiberGraphDistance || [];
    clear(ctx, canvas);

    if ($("showGhost").checked && frameIndex > 0) {
      for (const ids of fibers) {
        for (let j = 0; j < ids.length - 1; j++) line(ctx, xy(initial[ids[j]]), xy(initial[ids[j + 1]]), "#b9b1a5", .7, .45);
      }
    }

    fibers.forEach((ids, fid) => {
      const color = distances.length ? fiberColor(distances[fid]) : "#426f97";
      for (let j = 0; j < ids.length - 1; j++) line(ctx, xy(positions[ids[j]]), xy(positions[ids[j + 1]]), color, 1.1, .86);
    });

    if ($("showLinks").checked) {
      for (const linkDef of currentCase.crosslinks) {
        const a = xy(materialPoint(positions, linkDef, 0));
        const b = xy(materialPoint(positions, linkDef, 1));
        line(ctx, a, b, "#2a9c89", .75, .58);
        ctx.fillStyle = "#2a9c89";
        for (const p of [a, b]) {
          ctx.save(); ctx.translate(p[0], p[1]); ctx.rotate(Math.PI / 4); ctx.fillRect(-2, -2, 4, 4); ctx.restore();
        }
      }
    }

    if ($("showBeads").checked) {
      ctx.fillStyle = "#294f6b";
      for (let i = 0; i < positions.length; i++) {
        const p = xy(positions[i]);
        ctx.beginPath(); ctx.arc(p[0], p[1], 1.15, 0, Math.PI * 2); ctx.fill();
      }
    }
    ctx.fillStyle = "#17212b";
    for (const i of fixed) {
      const p = xy(positions[i]);
      ctx.fillRect(p[0] - 2, p[1] - 2, 4, 4);
    }

    const c = xy(center);
    const radiusValue = currentCase.mode === "cavity" ? metric.cavityRadius : metric.cellRadius;
    const radiusPx = (radiusValue || manifest.cellRadius) * (Math.min(canvas.width, canvas.height) - 48) / domain;
    ctx.fillStyle = currentCase.mode === "cavity" ? "rgba(128,85,163,.12)" : "rgba(212,94,56,.17)";
    ctx.strokeStyle = currentCase.mode === "cavity" ? "#8055a3" : "#a63e32";
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(c[0], c[1], radiusPx, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
    ctx.fillStyle = "#6e3d34";
    ctx.font = "600 13px system-ui";
    ctx.textAlign = "center";
    ctx.fillText(currentCase.mode === "cavity" ? "cavity" : "fixed-radius cell", c[0], c[1] + 4);

    if (caseId === "moving") {
      ctx.strokeStyle = "#a63e32"; ctx.lineWidth = 1.5; ctx.beginPath();
      currentCase.centers.slice(0, frameIndex + 1).forEach((point, index) => {
        const p = xy(point); if (index === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
      }); ctx.stroke();
    }

    const contacts = currentCase.contacts[frameIndex] || [];
    const total = currentCase.forces[frameIndex] || [];
    const parallel = currentCase.parallel[frameIndex] || [];
    const perpendicular = currentCase.perpendicular[frameIndex] || [];
    contacts.forEach((point, index) => {
      const p = xy(point);
      ctx.fillStyle = "#d45e38"; ctx.beginPath(); ctx.arc(p[0], p[1], 3.1, 0, 2 * Math.PI); ctx.fill();
      arrow(ctx, p, total[index] || [0, 0], "#a63e32", 2.0, 1.7);
      if ($("showComponents").checked) {
        arrow(ctx, p, parallel[index] || [0, 0], "#426f97", 1.2, 1.7);
        arrow(ctx, p, perpendicular[index] || [0, 0], "#8055a3", 1.2, 1.7);
      }
    });

    const amplification = Number($("amplification").value);
    ctx.strokeStyle = "rgba(38,48,55,.48)";
    for (let i = 0; i < positions.length; i += 18) {
      if (Math.hypot(initial[i][0], initial[i][1]) > manifest.cellRadius + 42) continue;
      const delta = [positions[i][0] - initial[i][0], positions[i][1] - initial[i][1]];
      if (Math.hypot(delta[0], delta[1]) < 1e-5) continue;
      const start = xy(initial[i]);
      const amplifiedEnd = [initial[i][0] + amplification * delta[0], initial[i][1] + amplification * delta[1]];
      const end = xy(amplifiedEnd);
      line(ctx, start, end, "#4f565b", .8, .55);
    }
    updateMetrics(currentCase, metric);
    drawPlot(currentCase);
  }

  function formatTime(seconds) {
    if (seconds >= 3600) return `${(seconds / 3600).toFixed(seconds % 3600 ? 1 : 0)} h`;
    if (seconds >= 60) return `${(seconds / 60).toFixed(seconds % 60 ? 1 : 0)} min`;
    return `${seconds.toFixed(seconds < 10 ? 1 : 0)} s`;
  }

  const interpretations = {
    passive: "Negative control: a fixed steric circle with no traction produces no collagen motion or alignment. This rules out geometry insertion as the cause.",
    cavity: "Mechanics benchmark only: the purple contractile cavity prescribes displacement and tests fiber rotation plus bending-to-stretching recruitment. It is not a shrinking tumor cell.",
    normal: "Primary fixed-cell test: every directly loaded point lies inside the finite contact band; Gaussian weights distribute a fixed total force only among those points. Indirect fibers move solely through network connections.",
    dipole: "Secondary same-total-force control: traction is restricted to two opposite sectors. It asks whether anisotropic loading recruits a longer tensile path; it is not the primary biological assumption.",
    "clutch-independent": "Redundancy control: individual clutches bind and rupture independently. Many molecular rupture events can occur while at least one clutch keeps the attachment site connected.",
    clutch: "Primary kinetic test: site load is shared by the bound clutches. When one ruptures, force per remaining clutch rises; complete site failure makes traction zero, resets extension and permits tangent-guided reselection.",
    moving: "Diagnostic released-cell run: the ECM mechanics and shared-load sites are unchanged; only cell translation is released under the equal-and-opposite reaction. No global polarization value is supplied.",
  };

  function updateMetrics(currentCase, metric) {
    $("timeLabel").textContent = formatTime(currentCase.times[frameIndex]);
    $("objectLabel").textContent = currentCase.centralObject;
    const delta = metric.deltaRadialOrder || [];
    $("nearMetric").textContent = delta.length ? Number(delta[0]).toFixed(4) : "0.0000";
    $("parallelMetric").textContent = `${Number(metric.parallelMagnitude || 0).toFixed(3)} nN`;
    $("perpMetric").textContent = `${Number(metric.perpendicularMagnitude || 0).toFixed(3)} nN`;
    $("strainMetric").textContent = Number(metric.meanAbsAxialStrain || 0).toExponential(2);
    $("bendMetric").textContent = `${Number(metric.meanBendingChange || 0).toExponential(2)} µm`;
    const bound = (metric.boundCount || []).reduce((a, b) => a + b, 0);
    $("boundMetric").textContent = metric.boundCount && metric.boundCount.length ? `${bound} across ${metric.boundCount.length} sites` : "not active";
    $("eventMetric").textContent = `${metric.ruptures || 0} / ${metric.siteFailures || 0}`;
    $("motionMetric").textContent = `${Number(currentCase.netDisplacement || 0).toFixed(3)} / ${Number(currentCase.pathLength || 0).toFixed(3)} µm`;
    $("interpretation").textContent = interpretations[caseId];
    const multiplier = Number($("amplification").value);
    $("scaleNote").innerHTML = `<b>Geometry is true scale.</b> Gray displacement arrows are ${multiplier}×; force arrows use a separate capped display scale.`;
  }

  function drawPlot(currentCase) {
    clear(pctx, plot);
    const times = currentCase.times;
    const metrics = currentCase.metrics;
    const left = 62, right = 62, top = 24, bottom = 47;
    const w = plot.width - left - right, h = plot.height - top - bottom;
    const maxT = Math.max(...times, 1);
    const alignment = metrics.map(m => Number((m.deltaRadialOrder || [0])[0] || 0));
    const traction = metrics.map(m => Number(m.distributedForce || (m.siteForce || []).reduce((a, b) => a + b, 0) || 0));
    const maxA = Math.max(.02, ...alignment.map(v => Math.abs(v)));
    const maxF = Math.max(1, ...traction);
    const X = t => left + t / maxT * w;
    const YA = v => top + h / 2 - v / maxA * h * .46;
    const YF = v => top + h - v / maxF * h;
    pctx.strokeStyle = "#bdb5a8"; pctx.lineWidth = 1;
    pctx.strokeRect(left, top, w, h);
    line(pctx, [left, YA(0)], [left + w, YA(0)], "#cec6ba", 1, 1);
    const trace = (values, y, color) => {
      pctx.strokeStyle = color; pctx.lineWidth = 2; pctx.beginPath();
      values.forEach((value, i) => { const px = X(times[i]), py = y(value); if (i === 0) pctx.moveTo(px, py); else pctx.lineTo(px, py); });
      pctx.stroke();
    };
    trace(alignment, YA, "#d45e38");
    trace(traction, YF, "#a63e32");
    pctx.fillStyle = "#26343c"; pctx.font = "13px system-ui";
    pctx.textAlign = "left"; pctx.fillText(`ΔSᵣ ±${maxA.toFixed(3)}`, 7, 18);
    pctx.fillStyle = "#d45e38"; pctx.fillText("near-cell ΔSᵣ", left + 8, top + 17);
    pctx.fillStyle = "#a63e32"; pctx.fillText("traction (nN)", left + 138, top + 17);
    pctx.fillStyle = "#26343c"; pctx.textAlign = "right"; pctx.fillText(`${maxF.toFixed(1)} nN`, plot.width - 5, 18);
    pctx.textAlign = "center"; pctx.fillText(formatTime(0), left, plot.height - 16); pctx.fillText(formatTime(maxT), left + w, plot.height - 16);
    const x = X(times[frameIndex]);
    line(pctx, [x, top], [x, top + h], "#17212b", 1, .7);
  }

  async function setCase(id) {
    $("interpretation").textContent = "Loading the selected precomputed trajectory…";
    try {
      await loadCase(id);
    } catch (error) {
      $("interpretation").textContent = `${error.message}. Open this page from GitHub Pages or the supplied local server.`;
      return;
    }
    caseId = id;
    frameIndex = 0;
    playing = false;
    if (timer) clearInterval(timer);
    $("play").textContent = "Play";
    const currentCase = cases[id];
    $("time").max = Math.max(0, currentCase.times.length - 1);
    $("time").value = 0;
    const event = events[id];
    $("eventPanel").classList.toggle("visible", Boolean(event && event.times && event.times.length));
    if (event && event.times.length) {
      $("eventTime").max = event.times.length - 1;
      const firstFailure = (event.eventLog || []).find(item => item.kind === "site_failure");
      const failureFrame = firstFailure
        ? event.times.reduce((best, time, frame) => Math.abs(time - firstFailure.time) < Math.abs(event.times[best] - firstFailure.time) ? frame : best, 0)
        : 0;
      $("eventTime").value = failureFrame;
      drawEvent(event, failureFrame);
    }
    drawNetwork();
  }

  function drawEvent(event, index) {
    clear(ectx, eventCanvas);
    const positions = eventPositions(event)[index];
    const center = event.centers[index];
    // The event file contains only the local collagen neighborhood.  Give it a
    // true-scale cell-local viewport instead of reusing the 180-µm overview
    // extent; otherwise a real clutch failure is visually sub-pixel.
    const localSpan = 2 * (manifest.cellRadius + manifest.config.contactWidth + 3);
    const toEventXY = point => xyLocal(point, eventCanvas, center, localSpan);
    const map = new Map(event.localBeads.map((global, local) => [global, local]));
    for (const edge of edges) {
      if (!map.has(edge[0]) || !map.has(edge[1])) continue;
      line(ectx, toEventXY(positions[map.get(edge[0])]), toEventXY(positions[map.get(edge[1])]), "#426f97", 1.2, .85);
    }
    positions.forEach(point => {
      const p = toEventXY(point); ectx.fillStyle = "#294f6b"; ectx.beginPath(); ectx.arc(p[0], p[1], 1.5, 0, 2 * Math.PI); ectx.fill();
    });
    const c = toEventXY(center);
    const radiusPx = manifest.cellRadius * (Math.min(eventCanvas.width, eventCanvas.height) - 56) / localSpan;
    ectx.fillStyle = "rgba(212,94,56,.17)"; ectx.strokeStyle = "#a63e32"; ectx.lineWidth = 2;
    ectx.beginPath(); ectx.arc(c[0], c[1], radiusPx, 0, 2 * Math.PI); ectx.fill(); ectx.stroke();
    const contacts = event.contacts[index] || [];
    const bound = event.bound[index] || [];
    contacts.forEach((point, site) => {
      const p = toEventXY(point);
      ectx.fillStyle = "#d45e38"; ectx.beginPath(); ectx.arc(p[0], p[1], 3, 0, 2 * Math.PI); ectx.fill();
      const state = bound[site] || [];
      state.forEach((on, clutch) => {
        if (!on) return;
        const angle = 2 * Math.PI * clutch / Math.max(state.length, 1);
        line(ectx, p, [p[0] + 7 * Math.cos(angle), p[1] + 7 * Math.sin(angle)], "#8055a3", 1.2, .8);
      });
      const breaks = (event.ruptures[index] || [])[site] || [];
      if (breaks.some(Boolean)) {
        ectx.strokeStyle = "#a63e32"; ectx.lineWidth = 2;
        ectx.beginPath(); ectx.moveTo(p[0] - 5, p[1] - 5); ectx.lineTo(p[0] + 5, p[1] + 5); ectx.moveTo(p[0] + 5, p[1] - 5); ectx.lineTo(p[0] - 5, p[1] + 5); ectx.stroke();
      }
    });
    const scalePx = 5 * (Math.min(eventCanvas.width, eventCanvas.height) - 56) / localSpan;
    ectx.strokeStyle = "#18242b"; ectx.lineWidth = 2;
    ectx.beginPath(); ectx.moveTo(24, eventCanvas.height - 22); ectx.lineTo(24 + scalePx, eventCanvas.height - 22); ectx.stroke();
    ectx.fillStyle = "#18242b"; ectx.font = "12px system-ui";
    ectx.fillText("5 µm", 24, eventCanvas.height - 28);
    const totalBound = (event.boundCount[index] || []).reduce((a, b) => a + b, 0);
    const force = (event.siteForce[index] || []).reduce((a, b) => a + b, 0);
    const failures = (event.eventLog || []).filter(item => item.kind === "site_failure" && Math.abs(item.time - event.times[index]) < .051).length;
    $("eventTimeLabel").textContent = formatTime(event.times[index]);
    $("eventBound").textContent = String(totalBound);
    $("eventForce").textContent = `${force.toFixed(3)} nN`;
    $("eventCounts").textContent = `${event.ruptureCount[index]} / ${event.bindCount[index]}`;
    $("eventFailures").textContent = String(failures);
  }

  function renderGates() {
    const items = [
      ["Passive control", manifest.gates.passiveZero, "Zero-force geometry produces negligible motion and ΔSᵣ."],
      ["Constant cell radius", manifest.gates.biologicalCellRadiusConstant, "Every biological run preserves the declared radius."],
      ["Mechanics identities", manifest.gates.mechanicsImplementation, "Force decomposition is finite and ECM/cell reaction balances."],
      ["Clutch event sequence", manifest.gates.clutchRuptureObserved, "Load-dependent rupture and a timestamped complete site failure are both observed."],
      ["Event microscope", manifest.gates.completeSiteFailureObserved, "The high-frequency window is anchored to a complete site failure, not an arbitrary frame."],
      ["Alignment vs passive", manifest.gates.alignmentExcludesPassive, "95% seed-ensemble interval must exclude zero."],
      ["Indirect transmission", manifest.gates.oneHopExceedsUnconnected, "One-hop displacement interval must exceed unconnected control."],
    ];
    const pending = [
      ["Radius sensitivity", "The 8/10/12-µm helper is implemented; the full seed ensemble is not yet run."],
      ["Domain convergence", "The 180/270/360-µm constant-density comparison is implemented; the full comparison is not yet run."],
      ["Timestep sensitivity", "The full 0.05/0.025-s comparison has not yet been run."],
    ];
    $("gateGrid").innerHTML =
      items.map(([name, pass, text]) => `<article class="gate ${pass ? "pass" : "fail"}"><strong>${pass ? "PASS" : "NOT YET"} · ${name}</strong>${text}</article>`).join("") +
      pending.map(([name, text]) => `<article class="gate pending"><strong>NOT RUN · ${name}</strong>${text}</article>`).join("");
    const e = manifest.ensemble;
    $("ensembleText").innerHTML = `<b>Seed ensemble:</b> near-cell mean ΔS<sub>r</sub> = ${Number(e.near_delta_mean).toFixed(4)}, 95% CI [${Number(e.near_delta_ci95[0]).toFixed(4)}, ${Number(e.near_delta_ci95[1]).toFixed(4)}]. One-hop minus unconnected displacement = ${Number(e.one_hop_advantage_mean).toExponential(2)} µm, 95% CI [${Number(e.one_hop_advantage_ci95[0]).toExponential(2)}, ${Number(e.one_hop_advantage_ci95[1]).toExponential(2)}]. A failed gate is preserved as a result, not hidden by animation tuning.`;
  }

  $("caseSelect").addEventListener("change", event => setCase(event.target.value));
  $("amplification").addEventListener("change", drawNetwork);
  ["showBeads", "showLinks", "showGhost", "showComponents"].forEach(id => $(id).addEventListener("change", drawNetwork));
  $("time").addEventListener("input", event => { frameIndex = Number(event.target.value); drawNetwork(); });
  $("step").addEventListener("click", () => {
    const frames = cases[caseId].times.length;
    frameIndex = (frameIndex + 1) % frames; $("time").value = frameIndex; drawNetwork();
  });
  $("reset").addEventListener("click", () => { frameIndex = 0; $("time").value = 0; drawNetwork(); });
  $("play").addEventListener("click", () => {
    playing = !playing; $("play").textContent = playing ? "Pause" : "Play";
    if (timer) clearInterval(timer);
    if (playing) timer = setInterval(() => {
      const frames = cases[caseId].times.length;
      frameIndex = (frameIndex + 1) % frames; $("time").value = frameIndex; drawNetwork();
    }, 450);
  });
  $("eventTime").addEventListener("input", event => drawEvent(events[caseId], Number(event.target.value)));

  const status = $("buildStatus");
  if (manifest.quick) {
    status.classList.add("quick");
    status.textContent = `QUICK PIPELINE CHECK · fixed ${formatTime(manifest.actualClocks.fixedSeconds)} · moving ${formatTime(manifest.actualClocks.movingSeconds)} · not the declared 2 h / 6 h evidence run`;
  } else {
    status.textContent = `FULL EVIDENCE BUILD · fixed ${formatTime(manifest.actualClocks.fixedSeconds)} · moving ${formatTime(manifest.actualClocks.movingSeconds)}`;
  }
  renderGates();
  setCase(caseId);
})();
