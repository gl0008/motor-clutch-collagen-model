(function () {
  'use strict';
  const D = window.G5_DATA;
  const $ = q => document.querySelector(q);
  if (!D) { $('#loadState').textContent = 'G5 data missing (run build_web.py).'; return; }

  const state = { mode: 'contract', frame: 0, playing: true, last: 0, view: 'full', links: false };

  // ---- decode helpers -----------------------------------------------------------
  function decodeFrames(pack, initialFlat) {
    const raw = atob(pack.b64), bytes = new Uint8Array(raw.length);
    for (let i = 0; i < raw.length; i++) bytes[i] = raw.charCodeAt(i);
    const q = new Int16Array(bytes.buffer), F = pack.F, N = pack.N, s = D.scale;
    const out = new Array(F);
    for (let f = 0; f < F; f++) {
      const row = new Float32Array(N * 2);
      for (let i = 0; i < N * 2; i++) row[i] = initialFlat[i] + q[f * N * 2 + i] * s;
      out[f] = row;
    }
    return out;
  }
  function prep(ds) {
    if (ds._pos) return ds;
    ds._pos = decodeFrames(ds.pos, ds.initial);
    ds._edges = ds.edges;
    // invade stores cells as an encoded pack (+cells0 base); contract stores a plain
    // fixed-cell array — only decode the pack.
    if (ds.cells && ds.cells.b64) ds._cells = decodeFrames(ds.cells, ds.cells0);
    return ds;
  }

  // ---- colour by radial order (blue tangential -> red radial) -------------------
  function coolwarm(t) { // t in [-1,1]
    t = (t + 1) / 2; t = Math.max(0, Math.min(1, t));
    const r = Math.round(59 + t * (180 - 59) + (t > .5 ? (t - .5) * 130 : 0));
    const g = Math.round(76 + (1 - Math.abs(t - .5) * 2) * 90);
    const b = Math.round(192 - t * (192 - 38));
    return `rgb(${Math.min(255, r)},${g},${Math.min(255, b)})`;
  }

  const cv = $('#stage'), ctx = cv.getContext('2d');

  function draw() {
    const ds = prep(D[state.mode]);
    const pos = ds._pos[Math.min(state.frame, ds._pos.length - 1)];
    const edges = ds._edges, R = D.cellRadius;
    const cells = state.mode === 'invade' ? ds._cells[Math.min(state.frame, ds._cells.length - 1)] : null;
    const cellsFlat = cells || ds.cells;         // Float32Array or plain array
    const nCells = cellsFlat.length / 2;
    const c0 = state.mode === 'invade' ? ds.cells0 : null;

    // world -> screen
    const span = state.view === 'full' ? D.domain : Math.min(D.domain, 190);
    const s = Math.min(cv.width, cv.height) / span;
    const T = (x, y) => [cv.width / 2 + x * s, cv.height / 2 - y * s];

    ctx.fillStyle = '#f4f0e7'; ctx.fillRect(0, 0, cv.width, cv.height);

    // fibres coloured by radial order of each segment
    for (let e = 0; e < edges.length; e += 2) {
      const i = edges[e], j = edges[e + 1];
      const ax = pos[i * 2], ay = pos[i * 2 + 1], bx = pos[j * 2], by = pos[j * 2 + 1];
      const dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy) || 1e-9;
      const mx = (ax + bx) / 2, my = (ay + by) / 2, mr = Math.hypot(mx, my) || 1e-9;
      const dot = (dx * mx + dy * my) / (L * mr);
      const order = 2 * dot * dot - 1;
      const A = T(ax, ay), B = T(bx, by);
      ctx.strokeStyle = coolwarm(order); ctx.lineWidth = 0.9;
      ctx.beginPath(); ctx.moveTo(A[0], A[1]); ctx.lineTo(B[0], B[1]); ctx.stroke();
    }

    // crosslinks (optional)
    if (state.links && D[state.mode].crosslinks) { /* reserved */ }

    // cell trails (invade)
    if (state.mode === 'invade' && c0) {
      ctx.strokeStyle = 'rgba(20,20,20,.35)'; ctx.lineWidth = 0.7;
      for (let k = 0; k < nCells; k++) {
        const a = T(c0[k * 2], c0[k * 2 + 1]), b = T(cellsFlat[k * 2], cellsFlat[k * 2 + 1]);
        ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
      }
    }

    // cell disks
    const rp = R * s;
    for (let k = 0; k < nCells; k++) {
      const p = T(cellsFlat[k * 2], cellsFlat[k * 2 + 1]);
      ctx.fillStyle = state.mode === 'invade' ? 'rgba(44,127,184,.88)' : 'rgba(80,80,80,.5)';
      ctx.beginPath(); ctx.arc(p[0], p[1], rp, 0, Math.PI * 2); ctx.fill();
    }

    // readouts
    const ds2 = D[state.mode];
    $('#tlabel').textContent = 't = ' + fmt(ds2.times[Math.min(state.frame, ds2.times.length - 1)]);
    if (state.mode === 'contract')
      $('#metric').innerHTML = 'network radial order <b>' + ds2.order[Math.min(state.frame, ds2.order.length - 1)].toFixed(3) + '</b> <span>(−1 tangential → +1 radial)</span>';
    else
      $('#metric').innerHTML = 'mean cell invasion <b>+' + ds2.invaded[Math.min(state.frame, ds2.invaded.length - 1)].toFixed(2) + ' µm</b> <span>(outward)</span>';
    $('#scrub').max = ds2.pos.F - 1;
    $('#scrub').value = state.frame;
  }
  function fmt(t) { return t >= 60 ? (t / 60).toFixed(t % 60 ? 1 : 0) + ' min' : t.toFixed(0) + ' s'; }

  function setMode(m) {
    state.mode = m; state.frame = 0;
    document.querySelectorAll('.modebtn').forEach(x => x.classList.toggle('active', x.dataset.mode === m));
    $('#caption').textContent = m === 'contract'
      ? 'Fixed organoid contracts; the connected near-field collagen reorganises toward radial.'
      : 'Released cells feel the reaction of their grip-and-reel traction and invade outward (trails).';
    draw();
  }

  document.querySelectorAll('.modebtn').forEach(b => b.onclick = () => setMode(b.dataset.mode));
  $('#play').onclick = () => { state.playing = !state.playing; $('#play').textContent = state.playing ? 'Pause' : 'Play'; };
  $('#view').onchange = e => { state.view = e.target.value; draw(); };
  $('#scrub').oninput = e => { state.playing = false; $('#play').textContent = 'Play'; state.frame = +e.target.value; draw(); };

  function tick(now) {
    if (state.playing && now - state.last > 120) {
      const F = D[state.mode].pos.F;
      state.frame = state.frame >= F - 1 ? 0 : state.frame + 1;
      draw(); state.last = now;
    }
    requestAnimationFrame(tick);
  }
  $('#loadState').textContent = 'loaded';
  setMode('contract');
  requestAnimationFrame(tick);
})();
