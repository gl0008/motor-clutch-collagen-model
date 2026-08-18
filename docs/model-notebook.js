(function(){
  const order=['g1-v0','g1-v1','g1-v2','g1-v3','g1-v4','g2-v2','g2-v3','g2-v4','g3-a','g3-b','g3-c'];
  const query=new URLSearchParams(location.search);
  const id=NOTEBOOK_MODELS[query.get('model')]?query.get('model'):'g2-v2';
  const model=NOTEBOOK_MODELS[id];
  const modelLabel=m=>m.generation==='G3'?m.version:`${m.generation} ${m.version}`;
  const esc=value=>String(value).replace(/[&<>"]/g,ch=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[ch]));
  const list=items=>`<ul>${items.map(x=>`<li>${x}</li>`).join('')}</ul>`;
  const button=(href,label,secondary=false)=>`<a class="button${secondary?' secondary':''}" href="${href}" target="_blank" rel="noopener">${label}</a>`;

  document.title=`${modelLabel(model)} · ${model.title}`;
  document.querySelector('#hero').classList.toggle('g2',model.generation==='G2');
  document.querySelector('#hero').classList.toggle('g3',model.generation==='G3');
  document.querySelector('#eyebrow').textContent=`${modelLabel(model)} · ${model.status}`;
  document.querySelector('#title').textContent=model.title;
  document.querySelector('#subtitle').textContent=model.subtitle;
  document.querySelector('#heroMeta').innerHTML=`<b>Question this version asks</b>${model.question}<br><br><b>Units</b>${model.units}`;
  document.querySelector('#versionNav').innerHTML=order.map(key=>{
    const m=NOTEBOOK_MODELS[key];
    return `<a class="${key===id?'current':''}" href="model-notebook.html?model=${key}">${modelLabel(m)}</a>`;
  }).join('');

  const brief=model.brief.map(x=>`<li>${x}</li>`).join('');
  const chain=model.chain.map((x,i)=>`${i?'<i>→</i>':''}<span>${x}</span>`).join('');
  let html=`
  <section class="section">
    <h2>Purpose of this version</h2>
    <p class="lead">${model.purpose}</p>
    <div class="brief ${model.generation==='G2'?'g2-accent':model.generation==='G3'?'g3-accent':''}"><div class="label">Your brief · condensed from our conversation</div><ul>${brief}</ul></div>
    <h3>Causal chain you wanted to isolate</h3><div class="causal-chain">${chain}</div>
    ${model.generation==='G3'?'<div class="professor-brief-link"><b>Presenting G3 tomorrow?</b><span>Open the one-page storyline with the model chain, key equations, paper map, current numbers and wording boundaries.</span><a class="button" href="g3-professor-brief.html">Open G3 professor briefing</a></div>':''}
  </section>`;

  if(model.simulation){
    html+=`<section class="section"><h2>Simulation</h2><p>${model.simulation.note}</p><div class="simulation-shell"><iframe title="${esc(model.title)} interactive simulation" src="${model.simulation.src}" loading="lazy"></iframe></div><div class="simulation-links">${button(model.simulation.src,'Open simulation full screen')}${model.simulation.gif?button(model.simulation.gif,'Open updated GIF',true):''}${button(model.source,'View model source',true)}</div></section>`;
  }else{
    html+=`<section class="section"><h2>Archived prototype artifacts</h2><div class="empty-sim"><h3>No standalone browser animation was published for V0</h3><p>The equations and code are preserved below. This is intentional: V0 is the question-forming prototype, not the accepted collagen-network baseline.</p></div><div class="archive-grid" style="margin-top:12px">${model.artifacts.map(x=>`<article class="archive-card"><b>${x.title}</b><p>${x.note}</p>${button(x.href,'Open source',true)}</article>`).join('')}</div></section>`;
  }

  const collectEquations=(key,trail=[])=>{
    const m=NOTEBOOK_MODELS[key];
    let groups=[];
    if(m.inherits)groups=groups.concat(collectEquations(m.inherits,[...trail,key]));
    return groups.concat((m.equationGroups||[]).map(group=>({...group,source:key===id?null:modelLabel(m)})));
  };
  const groups=collectEquations(id);
  html+=`<section class="section"><h2>Complete implemented equation set</h2><div class="scope-note"><b>Scope:</b> every governing, constitutive, stochastic, coupling, time-update and reported-metric equation that determines this version is shown. Low-level line-intersection and random-hash arithmetic are implementation utilities, not additional physics. Highlighted cards are the key equations to emphasize in a meeting.</div>`;
  html+=groups.map(group=>`<div class="equation-group"><div class="equation-group-head"><h3>${group.title}</h3>${group.source?`<span class="inherited">inherited from ${group.source}</span>`:''}</div>${group.intro?`<p>${group.intro}</p>`:''}<div class="equation-list">${group.equations.map(eq=>`<article class="equation-card ${eq.key?'key':''}">${eq.key?'<div class="equation-key-label">KEY EQUATION</div>':''}<h4>${eq.name}</h4><div class="math">${esc(eq.tex)}</div><p>${eq.explain}</p>${eq.units?`<small><b>Units/meaning:</b> ${eq.units}</small>`:''}</article>`).join('')}</div></div>`).join('');
  html+='</section>';

  html+=`<section class="section"><h2>How the simulation runs</h2><ol class="algorithm">${model.algorithm.map(x=>`<li>${x}</li>`).join('')}</ol></section>`;
  html+=`<section class="section"><h2>Paper support and evidence boundary</h2><div class="paper-grid">${model.papers.map(p=>`<article class="paper"><a href="${p.href}" target="_blank" rel="noopener">${p.title}</a><p>${p.support}</p><div class="boundary"><b>Does not prove:</b> ${p.boundary}</div></article>`).join('')}</div></section>`;
  html+=`<section class="section"><h2>What this version produced / should be checked</h2><div class="result-grid">${model.results.map(r=>`<article class="result-card"><b>${r.name}</b><p>${r.value}</p></article>`).join('')}</div></section>`;
  html+=`<section class="section ${model.generation==='G1'?'however':'limitation'}"><h2>${model.limitation.heading}</h2><p class="lead">${model.limitation.text}</p>${model.limitation.items?list(model.limitation.items):''}<p><b>Leads to:</b> ${model.limitation.leadsTo}</p></section>`;
  html+=`<section class="section talk"><h2>How to explain it to your professor</h2><p class="talk-script">${model.talk.script}</p><h3>If asked “what is the assumption?”</h3>${list(model.talk.assumptions)}<h3>If asked “what would falsify or improve this stage?”</h3>${list(model.talk.nextChecks)}</section>`;

  const index=order.indexOf(id),prev=order[index-1],next=order[index+1];
  html+=`<nav class="footer-nav">${prev?`<a href="model-notebook.html?model=${prev}">← ${modelLabel(NOTEBOOK_MODELS[prev])}</a>`:'<span></span>'}${next?`<a href="model-notebook.html?model=${next}">${modelLabel(NOTEBOOK_MODELS[next])} →</a>`:'<span></span>'}</nav>`;
  document.querySelector('#notebook').innerHTML=html;
  if(window.MathJax?.typesetPromise)window.MathJax.typesetPromise();
})();
