(function(){
  const groups=window.EVOLUTION_GROUPS;
  const versions=window.EVOLUTION_VERSIONS;
  const byGroup=Object.fromEntries(groups.map(group=>[group.id,versions.filter(version=>version.group===group.id)]));
  const groupMap=document.querySelector('#generationMap');
  const story=document.querySelector('#generationStory');
  const versionEvolution=document.querySelector('#versionEvolution');
  const filters=document.querySelector('#generationFilters');
  const list=document.querySelector('#versionList');
  const detail=document.querySelector('#versionDetail');
  let activeGroup='all';
  let activeVersion='g5-ablation';

  const escapeHtml=value=>String(value).replace(/[&<>"']/g,char=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[char]));
  const feedbackHtml=feedback=>feedback.length
    ? feedback.map(item=>`<div class="feedback-item"><div><b>${escapeHtml(item.who)}</b><time>${escapeHtml(item.date)}</time></div><p>${escapeHtml(item.text)}</p><small><b>Model response:</b> ${escapeHtml(item.impact)}</small></div>`).join('')
    : '<p class="empty-feedback">No version-specific professor feedback has been recorded yet.</p>';

  function renderMap(){
    const selectedGroup=activeGroup==='all'?versions.find(version=>version.id===activeVersion).group:activeGroup;
    groupMap.innerHTML=groups.map((group,index)=>{
      const node=`<button class="generation-node group-${group.id}" type="button" data-group="${group.id}" aria-pressed="${selectedGroup===group.id}" title="${escapeHtml(group.phase)}"><span>${group.label}</span><b>${escapeHtml(group.title)}</b></button>`;
      if(index===groups.length-1)return node;
      return `${node}<div class="generation-edge"><span class="edge-line" aria-hidden="true">→</span><p><b>Why it changed</b>${escapeHtml(group.edge)}</p></div>`;
    }).join('');
  }

  function renderStory(groupId){
    const group=groups.find(item=>item.id===groupId) || groups[groups.length-1];
    story.innerHTML=`
      <div class="story-title"><span>${group.label}</span><div><small>Selected generation</small><h3>${escapeHtml(group.title)}</h3></div></div>
      <p class="story-phase">Internal path: ${escapeHtml(group.phase)}</p>
      <div class="story-columns">
        <div><b>Inherited</b><p>${escapeHtml(group.inherited)}</p></div>
        <div class="limitation"><b>Limitation</b><p>${escapeHtml(group.limitation)}</p></div>
        <div class="question"><b>New question</b><p>${escapeHtml(group.question)}</p></div>
      </div>
      <div class="factor-row">${group.factors.map(factor=>`<span>${escapeHtml(factor)}</span>`).join('')}</div>`;
  }

  function renderVersionEvolution(groupId){
    const chain=byGroup[groupId];
    versionEvolution.innerHTML=`
      <div class="version-evolution-head">
        <b>Evolution within ${escapeHtml(groups.find(group=>group.id===groupId).label)}</b>
        <span>Scientific reasoning path, not raw Git commit order. Each arrow states why the next version was needed; select a version for its full record below.</span>
      </div>
      <div class="version-map-scroll">
        <div class="version-map-track">
          ${chain.map((version,index)=>{
            const node=`<button class="version-map-node" type="button" data-version="${version.id}" aria-pressed="${version.id===activeVersion}"><span class="version-map-label">${escapeHtml(version.label)}</span><b>${escapeHtml(version.title)}</b><small>${escapeHtml(version.status)}</small></button>`;
            if(index===chain.length-1)return node;
            return `${node}<div class="version-map-edge" title="${escapeHtml(version.next)}"><span aria-hidden="true">→</span><p><b>Why ${escapeHtml(chain[index+1].label)}?</b>${escapeHtml(version.next)}</p></div>`;
          }).join('')}
        </div>
      </div>`;
  }

  function renderFilters(){
    filters.innerHTML=[{id:'all',label:'All versions'},...groups].map(group=>`<button type="button" data-filter="${group.id}" aria-pressed="${activeGroup===group.id}">${group.label}</button>`).join('');
  }

  function renderList(){
    const visible=activeGroup==='all'?versions:byGroup[activeGroup];
    list.innerHTML=visible.map(version=>`
      <article class="version-row ${version.id===activeVersion?'selected':''}" data-row="${version.id}">
        <button class="version-select" type="button" data-version="${version.id}" aria-pressed="${version.id===activeVersion}">
          <span class="version-label">${escapeHtml(version.label)}</span>
          <span class="version-copy"><b>${escapeHtml(version.title)}</b><small>${escapeHtml(version.status)}</small></span>
          <span class="version-arrow">→</span>
        </button>
        <details class="feedback-fold">
          <summary>Professor feedback <span>${version.feedback.length}</span></summary>
          <div class="feedback-body">${feedbackHtml(version.feedback)}</div>
        </details>
      </article>`).join('');
  }

  function resolveQuestion(version){
    const notebook=window.NOTEBOOK_MODELS?.[version.id];
    return notebook?.question || version.change;
  }

  function renderDetail(versionId){
    const version=versions.find(item=>item.id===versionId) || versions[0];
    activeVersion=version.id;
    const notebook=window.NOTEBOOK_MODELS?.[version.id];
    const primaryHref=version.notebook || version.lab || version.source;
    const primaryLabel=version.notebook?'Open complete model notebook':version.lab?'Open simulation and full stage notes':'Open model source notes';
    const secondary=version.lab && version.notebook?`<a class="secondary-link" href="${version.lab}">Open interactive simulation</a>`:'';
    detail.innerHTML=`
      <div class="detail-head">
        <div><span class="detail-label">${escapeHtml(version.label)} · ${escapeHtml(version.status)}</span><h3>${escapeHtml(version.title)}</h3></div>
        <span class="record-id">${escapeHtml(version.id)}</span>
      </div>
      <div class="reasoning-band">
        <div><small>Previous limitation</small><p>${escapeHtml(version.limitation)}</p></div>
        <span>→</span>
        <div><small>Question this version asks</small><p>${escapeHtml(resolveQuestion(version))}</p></div>
      </div>
      <section class="detail-section"><h4>What changed / isolated</h4><p>${escapeHtml(version.change)}</p><div class="factor-row">${version.factors.map(factor=>`<span>${escapeHtml(factor)}</span>`).join('')}</div></section>
      <section class="detail-section equation-section"><h4>Key equation</h4><div class="equation">${version.equation}</div></section>
      <section class="detail-section provenance"><h4>Parameter provenance</h4><p>${escapeHtml(version.parameter)}</p></section>
      <section class="next-step"><small>Why the next version was needed</small><p>${escapeHtml(version.next)}</p></section>
      <details class="detail-feedback">
        <summary>Professor feedback <span>${version.feedback.length}</span></summary>
        <div class="feedback-body">${feedbackHtml(version.feedback)}</div>
      </details>
      <div class="detail-actions"><a class="primary-link" href="${primaryHref}">${primaryLabel}</a>${secondary}</div>
      ${notebook?`<p class="source-note">The complete page includes inherited equations, algorithm order, paper support, evidence boundaries, results and assumptions.</p>`:''}`;
    if(window.MathJax?.typesetPromise)window.MathJax.typesetPromise([detail]);
  }

  function selectGroup(groupId,focusFirst=true){
    activeGroup=groupId;
    if(groupId!=='all' && (focusFirst || !byGroup[groupId].some(version=>version.id===activeVersion))){
      activeVersion=byGroup[groupId][0].id;
    }
    const selectedGroup=groupId==='all'?versions.find(version=>version.id===activeVersion).group:groupId;
    renderMap();
    renderStory(selectedGroup);
    renderVersionEvolution(selectedGroup);
    renderFilters();
    renderList();
    renderDetail(activeVersion);
  }

  document.addEventListener('click',event=>{
    const groupButton=event.target.closest('[data-group]');
    if(groupButton){selectGroup(groupButton.dataset.group);return;}
    const filterButton=event.target.closest('[data-filter]');
    if(filterButton){selectGroup(filterButton.dataset.filter,false);return;}
    const versionButton=event.target.closest('[data-version]');
    if(versionButton){
      activeVersion=versionButton.dataset.version;
      document.querySelectorAll('.version-row').forEach(row=>row.classList.toggle('selected',row.dataset.row===activeVersion));
      document.querySelectorAll('[data-version]').forEach(button=>button.setAttribute('aria-pressed',button.dataset.version===activeVersion));
      document.querySelectorAll('.version-map-node').forEach(button=>button.classList.toggle('selected',button.dataset.version===activeVersion));
      if(activeGroup==='all'){
        const selected=versions.find(version=>version.id===activeVersion);
        renderMap();
        renderStory(selected.group);
        renderVersionEvolution(selected.group);
      }
      renderDetail(activeVersion);
      if(versionButton.classList.contains('version-map-node')){
        detail.scrollIntoView({behavior:'smooth',block:'start'});
      }
      return;
    }
  });

  selectGroup('g5',false);
})();
