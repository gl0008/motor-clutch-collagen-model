(function(){
  const groups=window.EVOLUTION_GROUPS || [];
  const versions=window.EVOLUTION_VERSIONS || [];
  const logs=window.EXPERIMENT_LOG || [];
  const branches=window.REPOSITORY_BRANCHES || [];
  const evidenceByVersion=window.VERSION_EVIDENCE || {};
  const byGroup=Object.fromEntries(groups.map(group=>[group.id,versions.filter(version=>version.group===group.id)]));
  const groupMap=document.querySelector('#generationMap');
  const story=document.querySelector('#generationStory');
  const versionEvolution=document.querySelector('#versionEvolution');
  const filters=document.querySelector('#generationFilters');
  const list=document.querySelector('#versionList');
  const detail=document.querySelector('#versionDetail');
  const branchRegistry=document.querySelector('#branchRegistry');
  const experimentLog=document.querySelector('#experimentLog');
  let activeMapGroup='g5';
  let activeFilter='all';
  let activeVersion='g5-ablation';

  const repository='https://github.com/gl0008/motor-clutch-collagen-model';
  const escapeHtml=value=>String(value ?? '').replace(/[&<>"']/g,char=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[char]));
  const branchHref=branch=>`${repository}/tree/${branch}`;
  const commitHref=commit=>`${repository}/commit/${commit}`;
  const feedbackHtml=feedback=>feedback.length
    ? feedback.map(item=>`<div class="feedback-item"><div><b>${escapeHtml(item.who)}</b><time>${escapeHtml(item.date)}</time></div><p>${escapeHtml(item.text)}</p><small><b>Model response:</b> ${escapeHtml(item.impact)}</small></div>`).join('')
    : '<p class="empty-feedback">No version-specific professor feedback has been recorded yet.</p>';

  function renderMap(){
    groupMap.innerHTML=groups.map((group,index)=>{
      const node=`<button class="generation-node group-${group.id}" type="button" data-group="${group.id}" aria-pressed="${activeMapGroup===group.id}" title="${escapeHtml(group.phase)}"><span>${group.label}</span><b>${escapeHtml(group.title)}</b></button>`;
      if(index===groups.length-1)return node;
      return `${node}<div class="generation-edge"><span class="edge-line" aria-hidden="true">→</span><p><b>Limitation → next question</b>${escapeHtml(group.edge)}</p></div>`;
    }).join('');
  }

  function renderStory(groupId){
    const group=groups.find(item=>item.id===groupId) || groups[0];
    if(!group)return;
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
    const group=groups.find(item=>item.id===groupId);
    const chain=byGroup[groupId] || [];
    if(!group)return;
    versionEvolution.innerHTML=`
      <div class="version-evolution-head">
        <b>Evolution within ${escapeHtml(group.label)}</b>
        <span>This follows the scientific reasoning, not raw Git order. Select a version to open its evidence-linked record.</span>
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
    filters.innerHTML=[{id:'all',label:'All versions'},...groups].map(group=>`<button type="button" data-filter="${group.id}" aria-pressed="${activeFilter===group.id}">${escapeHtml(group.label)}</button>`).join('');
  }

  function renderList(){
    const visible=activeFilter==='all'?versions:(byGroup[activeFilter] || []);
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

  function theoryLink(groupId){
    const files={g3:'g3_model_theory.pdf',g4:'g4_model_theory.pdf',g5:'g5_model_theory.pdf'};
    return files[groupId]?`<a class="secondary-link" href="theory/${files[groupId]}">Open ${groupId.toUpperCase()} theory guide</a>`:'';
  }

  function renderDetail(versionId){
    const version=versions.find(item=>item.id===versionId) || versions[0];
    if(!version)return;
    activeVersion=version.id;
    const notebook=window.NOTEBOOK_MODELS?.[version.id];
    const evidence=evidenceByVersion[version.id];
    const primaryHref=version.notebook || version.lab || version.source;
    const primaryLabel=version.notebook?'Open complete model notebook':version.lab?'Open simulation and stage notes':'Open model source notes';
    const simulation=version.lab && version.notebook?`<a class="secondary-link" href="${version.lab}">Open simulation</a>`:'';
    const primary=primaryHref?`<a class="primary-link" href="${primaryHref}">${primaryLabel}</a>`:'';
    const repositoryEvidence=evidence?`
      <section class="detail-section repository-evidence">
        <h4>Repository evidence</h4>
        <p>This placement is anchored to the preserved branch and implementation commit below.</p>
        <div class="evidence-links">
          <a href="${branchHref(evidence.branch)}"><span>Branch</span><code>${escapeHtml(evidence.branch)}</code></a>
          <a href="${commitHref(evidence.commit)}"><span>Commit</span><code>${escapeHtml(evidence.commit)}</code></a>
          <a href="#experiments" data-view-link="experiments"><span>Human record</span><b>Find the experiment logic</b></a>
        </div>
      </section>`:'';
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
      <section class="detail-section"><h4>What changed or was isolated</h4><p>${escapeHtml(version.change)}</p><div class="factor-row">${version.factors.map(factor=>`<span>${escapeHtml(factor)}</span>`).join('')}</div></section>
      <section class="detail-section equation-section"><h4>Key equation</h4><div class="equation">${version.equation}</div></section>
      <section class="detail-section provenance"><h4>Parameter provenance</h4><p>${escapeHtml(version.parameter)}</p></section>
      ${repositoryEvidence}
      <section class="next-step"><small>Why the next version was needed</small><p>${escapeHtml(version.next)}</p></section>
      <details class="detail-feedback">
        <summary>Professor feedback <span>${version.feedback.length}</span></summary>
        <div class="feedback-body">${feedbackHtml(version.feedback)}</div>
      </details>
      <div class="detail-actions">${primary}${simulation}${theoryLink(version.group)}</div>
      ${notebook?'<p class="source-note">The complete notebook includes inherited equations, algorithm order, paper support, evidence boundaries, results and assumptions.</p>':''}`;
    if(window.MathJax?.typesetPromise)window.MathJax.typesetPromise([detail]);
  }

  function renderBranches(){
    branchRegistry.innerHTML=branches.map(branch=>`
      <article class="branch-card">
        <div><span class="branch-state state-${escapeHtml(branch.state).replace(/\s+/g,'-')}">${escapeHtml(branch.state)}</span><b>${escapeHtml(branch.role)}</b></div>
        <a href="${branchHref(branch.name)}"><code>${escapeHtml(branch.name)}</code></a>
        <p>${escapeHtml(branch.description)}</p>
      </article>`).join('');
  }

  function renderExperiments(){
    experimentLog.innerHTML=logs.map((entry,index)=>`
      <article class="log-entry">
        <header class="log-head">
          <div><span class="log-index">${String(logs.length-index).padStart(2,'0')}</span><div><small>${escapeHtml(entry.date)} · ${escapeHtml(entry.generation)}</small><h3>${escapeHtml(entry.version)}</h3></div></div>
          <div class="log-repository"><a href="${branchHref(entry.branch)}"><code>${escapeHtml(entry.branch)}</code></a><a href="${commitHref(entry.commit)}"><code>${escapeHtml(entry.commit)}</code></a></div>
        </header>
        <p class="log-parent"><b>Inherited from:</b> ${escapeHtml(entry.parent)}</p>
        <div class="log-flow">
          <div><small>Question</small><p>${escapeHtml(entry.question)}</p></div>
          <span aria-hidden="true">→</span>
          <div><small>Change</small><p>${escapeHtml(entry.change)}</p></div>
          <span aria-hidden="true">→</span>
          <div><small>Result</small><p>${escapeHtml(entry.result)}</p></div>
        </div>
        <div class="log-followup"><div><small>Limitation</small><p>${escapeHtml(entry.limitation)}</p></div><span aria-hidden="true">→</span><div><small>Next question</small><p>${escapeHtml(entry.next)}</p></div></div>
        <div class="evidence-links compact"><span>Supporting evidence</span>${entry.evidence.map(item=>`<a href="${item.href}">${escapeHtml(item.label)}</a>`).join('')}</div>
      </article>`).join('');
  }

  function selectGroup(groupId){
    if(!byGroup[groupId])return;
    activeMapGroup=groupId;
    if(!byGroup[groupId].some(version=>version.id===activeVersion))activeVersion=byGroup[groupId][0].id;
    renderMap();
    renderStory(groupId);
    renderVersionEvolution(groupId);
  }

  function selectFilter(groupId){
    activeFilter=groupId;
    const visible=groupId==='all'?versions:(byGroup[groupId] || []);
    if(!visible.some(version=>version.id===activeVersion))activeVersion=visible[0]?.id || versions[0]?.id;
    renderFilters();
    renderList();
    renderDetail(activeVersion);
  }

  function selectVersion(versionId,openLibrary=false){
    const version=versions.find(item=>item.id===versionId);
    if(!version)return;
    activeVersion=version.id;
    if(openLibrary){
      activeFilter=version.group;
      setView('versions',true);
    }
    renderFilters();
    renderList();
    renderDetail(version.id);
    renderVersionEvolution(activeMapGroup);
  }

  function setView(view,updateHash=false){
    const valid=['evolution','versions','experiments'];
    const selected=valid.includes(view)?view:'evolution';
    document.querySelectorAll('[data-view]').forEach(section=>{section.hidden=section.dataset.view!==selected;});
    document.querySelectorAll('[data-view-link]').forEach(link=>{
      const current=link.dataset.viewLink===selected;
      link.classList.toggle('current',current);
      if(current)link.setAttribute('aria-current','page'); else link.removeAttribute('aria-current');
    });
    if(updateHash && location.hash!==`#${selected}`)location.hash=selected;
  }

  document.addEventListener('click',event=>{
    const viewLink=event.target.closest('[data-view-link]');
    if(viewLink){event.preventDefault();setView(viewLink.dataset.viewLink,true);return;}
    const groupButton=event.target.closest('[data-group]');
    if(groupButton){selectGroup(groupButton.dataset.group);return;}
    const filterButton=event.target.closest('[data-filter]');
    if(filterButton){selectFilter(filterButton.dataset.filter);return;}
    const versionButton=event.target.closest('[data-version]');
    if(versionButton){selectVersion(versionButton.dataset.version,versionButton.classList.contains('version-map-node'));return;}
  });

  window.addEventListener('hashchange',()=>setView(location.hash.slice(1)));
  renderMap();
  renderStory(activeMapGroup);
  renderVersionEvolution(activeMapGroup);
  renderFilters();
  renderList();
  renderDetail(activeVersion);
  renderBranches();
  renderExperiments();
  setView(location.hash.slice(1));
})();
