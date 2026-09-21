/* G5 organoid interactive viewer — mirrors gl0008 g4-v4e.js (canvas 2D, packed int16 deltas).
   Loads window.G5_MANIFEST + on-demand window.G5_CASES[tag]; redraws every frame at true scale. */
(function () {
  "use strict";
  var $ = function (id) { return document.getElementById(id); };
  var manifest = window.G5_MANIFEST || { cases: [] };
  var C = null;           // current case (raw)
  var POS = null;         // decoded positions: Array(nFrames) of Float32Array(2*nBeads)
  var EDGES = null;       // Int32Array pairs
  var frame = 0, playing = false, lastT = 0, acc = 0;

  var COL = {
    bg: "#f7f4ec", fibre: "#476f91", ghost: "#b7afa2", link: "#2b9b83",
    trail: "#9a8f80", arrow: "#17232b",
    cohFill: "rgba(127,168,201,.55)", cohEdge: "#2c5f80",
    emtFill: "rgba(217,138,122,.65)", emtEdge: "#a83d35",
  };

  function unpack(pack, base) {
    var bin = atob(pack.base64), bytes = new Uint8Array(bin.length), i;
    for (i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    var v = new Int16Array(bytes.buffer), nf = pack.shape[0], np = pack.shape[1], sc = pack.scale;
    var out = new Array(nf), k = 0, f, j;
    for (f = 0; f < nf; f++) {
      var a = new Float32Array(2 * np);
      for (j = 0; j < np; j++) { a[2 * j] = base[j][0] + sc * v[k++]; a[2 * j + 1] = base[j][1] + sc * v[k++]; }
      out[f] = a;
    }
    return out;
  }

  function pxPerUm(canvas) { return (Math.min(canvas.width, canvas.height) - 40) / (2 * C.span); }
  function X(canvas, x) { return canvas.width / 2 + x * pxPerUm(canvas); }
  function Y(canvas, y) { return canvas.height / 2 - y * pxPerUm(canvas); }

  function redSet(f) {
    if (C.redPerFrame) return new Set(C.redPerFrame[f]);
    return new Set(C.redIds || []);
  }

  function draw() {
    var canvas = $("view"), ctx = canvas.getContext("2d"), f = frame;
    var s = pxPerUm(canvas), pos = POS[f], W = canvas.width, H = canvas.height;
    ctx.fillStyle = COL.bg; ctx.fillRect(0, 0, W, H);

    // initial ghost fibres
    if ($("ghost").checked && f > 0) {
      ctx.strokeStyle = COL.ghost; ctx.lineWidth = 0.5; ctx.globalAlpha = 0.5; ctx.beginPath();
      var b = C.base;
      for (var e = 0; e < EDGES.length; e += 2) {
        var i0 = EDGES[e], i1 = EDGES[e + 1];
        ctx.moveTo(X(canvas, b[i0][0]), Y(canvas, b[i0][1]));
        ctx.lineTo(X(canvas, b[i1][0]), Y(canvas, b[i1][1]));
      }
      ctx.stroke(); ctx.globalAlpha = 1;
    }

    // fibres (batched single path for performance)
    ctx.strokeStyle = COL.fibre; ctx.lineWidth = 0.7; ctx.globalAlpha = 0.72; ctx.beginPath();
    for (var g = 0; g < EDGES.length; g += 2) {
      var a0 = EDGES[g], a1 = EDGES[g + 1];
      ctx.moveTo(W / 2 + pos[2 * a0] * s, H / 2 - pos[2 * a0 + 1] * s);
      ctx.lineTo(W / 2 + pos[2 * a1] * s, H / 2 - pos[2 * a1 + 1] * s);
    }
    ctx.stroke(); ctx.globalAlpha = 1;

    // bead dots along fibres (Gloria's beaded-fibre texture; every 2nd bead)
    ctx.fillStyle = COL.fibre; ctx.globalAlpha = 0.75;
    for (var bd = 0; bd < pos.length; bd += 4) {
      ctx.fillRect(W / 2 + pos[bd] * s - 0.55, H / 2 - pos[bd + 1] * s - 0.55, 1.1, 1.1);
    }
    ctx.globalAlpha = 1;

    // crosslinks (dots on the moving fibres)
    if ($("links").checked) {
      ctx.fillStyle = COL.link;
      for (var x = 0; x < C.xlEdge.length; x++) {
        var ed = C.xlEdge[x], al = C.xlAlpha[x], e0 = EDGES[2 * ed], e1 = EDGES[2 * ed + 1];
        var px = pos[2 * e0] * (1 - al) + pos[2 * e1] * al, py = pos[2 * e0 + 1] * (1 - al) + pos[2 * e1 + 1] * al;
        var cx = W / 2 + px * s, cy = H / 2 - py * s;
        ctx.fillRect(cx - 1, cy - 1, 2, 2);
      }
    }

    // cell paths (trails t0..f)
    if ($("trails").checked && f > 0) {
      ctx.strokeStyle = COL.trail; ctx.lineWidth = 1; ctx.globalAlpha = 0.7;
      for (var ci = 0; ci < C.nCells; ci++) {
        ctx.beginPath();
        for (var t = 0; t <= f; t++) {
          var q = C.cells[t][ci], qx = X(canvas, q[0]), qy = Y(canvas, q[1]);
          if (t === 0) ctx.moveTo(qx, qy); else ctx.lineTo(qx, qy);
        }
        ctx.stroke();
      }
      ctx.globalAlpha = 1;
    }

    // net traction arrows (one per cell, fixed length, directional)
    if ($("arrows").checked) {
      var L = 1.5 * C.cellRadius * s;
      for (var k = 0; k < C.nCells; k++) {
        var nv = C.net[f][k], mag = Math.hypot(nv[0], nv[1]); if (mag < 1e-9) continue;
        var cc = C.cells[f][k], ox = X(canvas, cc[0]), oy = Y(canvas, cc[1]);
        var ex = ox + nv[0] / mag * L, ey = oy - nv[1] / mag * L;
        ctx.strokeStyle = COL.arrow; ctx.lineWidth = 1.6; ctx.beginPath();
        ctx.moveTo(ox, oy); ctx.lineTo(ex, ey); ctx.stroke();
        var ang = Math.atan2(ey - oy, ex - ox), hs = 5;
        ctx.fillStyle = COL.arrow; ctx.beginPath(); ctx.moveTo(ex, ey);
        ctx.lineTo(ex - hs * Math.cos(ang - .5), ey - hs * Math.sin(ang - .5));
        ctx.lineTo(ex - hs * Math.cos(ang + .5), ey - hs * Math.sin(ang + .5));
        ctx.fill();
      }
    }

    // cells
    var red = redSet(f), rp = C.cellRadius * s;
    for (var m = 0; m < C.nCells; m++) {
      var c = C.cells[f][m], isR = red.has(m), cx2 = X(canvas, c[0]), cy2 = Y(canvas, c[1]);
      ctx.fillStyle = isR ? COL.emtFill : COL.cohFill;
      ctx.strokeStyle = isR ? COL.emtEdge : COL.cohEdge; ctx.lineWidth = 1.6;
      ctx.beginPath(); ctx.arc(cx2, cy2, rp, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
    }

    // scale bar (50 um)
    var barPx = 50 * s, bx = 30, by = H - 28;
    ctx.strokeStyle = "#17232b"; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(bx, by); ctx.lineTo(bx + barPx, by); ctx.stroke();
    ctx.fillStyle = "#17232b"; ctx.font = "12px system-ui"; ctx.fillText("50 µm", bx + barPx / 2 - 16, by - 6);

    $("time").value = f; $("timeLabel").textContent = (C.times[f] / 60).toFixed(0) + " min";
    var inv = meanInv(f);
    $("metrics").innerHTML = "<b>" + C.label + "</b> — t = " + (C.times[f] / 60).toFixed(0)
      + " min &nbsp;·&nbsp; mean invasion " + (inv >= 0 ? "+" : "") + inv.toFixed(1) + " µm"
      + " &nbsp;·&nbsp; force-pair " + (C.summary.fpair || 0).toExponential(0)
      + (C.redPerFrame ? " &nbsp;·&nbsp; red = current leader (baton)" : (red.size ? " &nbsp;·&nbsp; red = partial-EMT/leader" : ""));
  }

  function meanInv(f) {
    var s = 0, c0 = C.cells[0], cf = C.cells[f];
    for (var i = 0; i < C.nCells; i++)
      s += Math.hypot(cf[i][0], cf[i][1]) - Math.hypot(c0[i][0], c0[i][1]);
    return s / C.nCells;
  }

  function tick(ts) {
    if (!playing) return;
    if (!lastT) lastT = ts;
    acc += ts - lastT; lastT = ts;
    var step = parseFloat($("speed").value);
    while (acc >= step) { acc -= step; frame++; if (frame >= C.nFrames) { frame = 0; } draw(); }
    requestAnimationFrame(tick);
  }

  function play() { playing = !playing; $("play").textContent = playing ? "Pause" : "Play"; lastT = 0; acc = 0; if (playing) requestAnimationFrame(tick); }

  function setup(tag) {
    C = window.G5_CASES[tag];
    POS = unpack(C.positions, C.base);
    EDGES = new Int32Array(C.edges.length * 2);
    for (var i = 0; i < C.edges.length; i++) { EDGES[2 * i] = C.edges[i][0]; EDGES[2 * i + 1] = C.edges[i][1]; }
    frame = 0; playing = false; $("play").textContent = "Play";
    $("time").max = C.nFrames - 1; $("time").value = 0;
    $("status").textContent = C.nFrames + " frames · " + C.nCells + " cells · " + C.base.length + " beads";
    draw();
  }

  function loadCase(tag) {
    $("status").textContent = "Loading " + tag + "…";
    if (window.G5_CASES && window.G5_CASES[tag]) { setup(tag); return; }
    var sc = document.createElement("script");
    sc.src = "g5-data/" + tag + ".js";
    sc.onload = function () { setup(tag); };
    sc.onerror = function () { $("status").textContent = "Failed to load " + tag; };
    document.body.appendChild(sc);
  }

  function init() {
    var sel = $("mode");
    (manifest.cases || []).forEach(function (m) {
      var o = document.createElement("option"); o.value = m.tag; o.textContent = m.label; sel.appendChild(o);
    });
    sel.onchange = function () { loadCase(sel.value); };
    $("play").onclick = play;
    $("step").onclick = function () { playing = false; $("play").textContent = "Play"; frame = (frame + 1) % C.nFrames; draw(); };
    $("reset").onclick = function () { playing = false; $("play").textContent = "Play"; frame = 0; draw(); };
    $("time").oninput = function () { playing = false; $("play").textContent = "Play"; frame = parseInt($("time").value, 10); draw(); };
    ["ghost", "trails", "links", "arrows"].forEach(function (id) { $(id).onchange = draw; });
    if ((manifest.cases || []).length) loadCase(manifest.cases[0].tag);
    else $("status").textContent = "No data — run export_web_viewer.py";
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
