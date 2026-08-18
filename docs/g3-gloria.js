(function(){
  const P=window.G3_WEB_DATA,NS='http://www.w3.org/2000/svg';
  const stageKey=new URLSearchParams(location.search).get('stage')||'g3a';
  const stages={
    g3a:{title:'G3A-R · protrusion-tip attachment and collagen loading',subtitle:'Both panels read the same solved 99-fibre trajectory. The left isolates growth and material-point capture; the right exposes the force-driven fibre response.',gif:'assets/g3_revision/g3a_multifibre_protrusion.gif',left:['Protrusion growth and attachment','g3a','mechanism'],right:['Force-driven fibre response','g3a','response'],view:'near',vectors:true,equation:'ζᵦ ṙᵢ = Fᵢˢᵗʳᵉᵗᶜʰ + Fᵢᵇᵉⁿᵈ + Fᵢˣˡⁱⁿᵏ + Fᵢᶜˡᵘᵗᶜʰ + Fᵢᶜᵒⁿᵗᵃᶜᵗ',what:'A green protrusion grows at finite speed. It can bind only when its explicit tip lies within 0.5 µm of a continuous collagen segment. Red clutch spokes then transmit tension to that stored material point; orange vectors and the recorded angle metric expose the resulting fibre response.'},
    g3b:{title:'G3B-R · cell-intrinsic polarity with adhesion feedback',subtitle:'The two panels use the same solved seed and aligned 99-fibre network. Only the post-attachment adhesion term is removed on the left.',gif:'assets/g3_revision/g3b_intrinsic_polarity_feedback.gif',left:['Adhesion feedback OFF','g3b_off','standard'],right:['Adhesion feedback ON','g3b_on','standard'],view:'near',vectors:false,equation:'aₛⁿ⁺¹ = 𝒩₊[aₛⁿ + (Δt/τₐ)(Acell πₛ − aₛⁿ) + ηₛ],   Σₛaₛ = 1',what:'The activity pool is conserved and cell intrinsic; no 0.65 direction or collagen angle is supplied. Bright green protrusions are active, pale green ones are retracting. In this matched short run, OFF and ON follow the same sector sequence, so feedback stabilization is not established.'},
    g3c:{title:'G3C-R · fixed versus reaction-driven released cell',subtitle:'Both panels use the same solved polarity, clutch and 99-fibre ECM loop. Only the right cell integrates the equal-and-opposite clutch/contact reaction.',gif:'assets/g3_revision/g3c_fixed_released.gif',left:['Fixed cell','g3c_fixed','standard'],right:['Released cell','g3c_released','standard'],view:'full',vectors:false,equation:'ṙc = Fcell/ζc,   φ̇c = τcell/ζr,   Fcell = −Σc∈bound fc',what:'The fixed cell records force but cannot translate or rotate. The released cell updates x, y and body angle from the solved reaction force and torque. There is no prescribed velocity or browser-side motion; the current displacement is real but very small.'}
  };
  const spec=stages[stageKey]||stages.g3a;
  let frame=0,playing=true,last=0,panels=[];

  function el(tag,attrs={}){const node=document.createElementNS(NS,tag);for(const[k,v]of Object.entries(attrs))node.setAttribute(k,v);return node}
  function decode(data){if(data.q)return;const enc=data.position_encoding,raw=atob(enc.base64),bytes=new Uint8Array(raw.length);for(let i=0;i<raw.length;i++)bytes[i]=raw.charCodeAt(i);data.q=new Int16Array(bytes.buffer)}
  function xy(data,f,i){const n=data.position_encoding.shape[1],k=(f*n+i)*2,q=data.position_encoding.quantum_um,p=data.initial_positions[i];return[p[0]+q*data.q[k],p[1]+q*data.q[k+1]]}
  function material(data,f,link,side){const edge=side==='a'?link.edge_a:link.edge_b,a=side==='a'?link.alpha_a:link.alpha_b,ij=data.edges[edge],p=xy(data,f,ij[0]),q=xy(data,f,ij[1]);return[(1-a)*p[0]+a*q[0],(1-a)*p[1]+a*q[1]]}
  function setup(svgId,metricId,entry,color){
    const [label,key,role]=entry,data=P[key],scene=document.querySelector('#'+svgId);decode(data);
    const defs=el('defs'),marker=el('marker',{id:'arrow-'+svgId,viewBox:'0 0 6 6',refX:5.25,refY:3,markerWidth:2.3,markerHeight:2.3,orient:'auto'});marker.append(el('path',{d:'M0,0 L6,3 L0,6 Z',fill:'#d45f32'}));defs.append(marker);scene.append(defs);
    const ghost=el('g'),bonds=el('g'),beads=el('g'),links=el('g'),vectors=el('g'),activity=el('g'),protrusions=el('g'),clutches=el('g'),anchors=el('g');scene.append(ghost,bonds,beads,links,vectors,activity,protrusions,clutches,anchors);
    data.fibers.forEach(ids=>ghost.append(el('polyline',{class:'ghost',points:ids.map(i=>data.initial_positions[i].join(',')).join(' ')})));
    const edgeEls=data.edges.map(()=>{const node=el('line',{class:'bond'});bonds.append(node);return node});
    const beadEls=data.initial_positions.map(()=>{const node=el('circle',{r:.30,class:'bead'});beads.append(node);return node});
    data.fixed.forEach((fixed,i)=>{if(fixed){const node=el('rect',{width:.75,height:.75,class:'anchor'});node.dataset.i=i;anchors.append(node)}});
    const initialCell=el('circle',{cx:data.cell_center[0][0],cy:data.cell_center[0][1],r:data.config.cell_radius,class:'initialCell'}),trail=el('polyline',{class:'trail'}),cell=el('circle',{r:data.config.cell_radius,class:'cell',stroke:color});scene.append(initialCell,trail,cell);
    const scale=el('g');scale.append(el('line',{class:'scale'}),el('text',{class:'scaleText'}));scene.append(scale);
    document.querySelector('#'+svgId+'Title').textContent=label;document.querySelector('#'+svgId+'Title').style.color=color;
    return{scene,metricId,data,role,color,edgeEls,beadEls,links,vectors,activity,protrusions,clutches,anchors,cell,initialCell,trail,scale};
  }
  function setView(panel,f){
    const data=panel.data,c=data.cell_center[f],full=document.querySelector('#view').value==='full';let box;
    if(full){const h=data.config.domain_size/2;box=[-h,-h,2*h,2*h]}else{const half=panel.role==='response'?38:28;box=[c[0]-half,c[1]-half,2*half,2*half]}
    panel.scene.setAttribute('viewBox',box.join(' '));const line=panel.scale.querySelector('line'),text=panel.scale.querySelector('text'),sx=box[0]+.055*box[2],sy=box[1]+.94*box[3];line.setAttribute('x1',sx);line.setAttribute('x2',sx+10);line.setAttribute('y1',sy);line.setAttribute('y2',sy);text.setAttribute('x',sx);text.setAttribute('y',sy-1.2);text.textContent='10 µm';
  }
  function metrics(panel,f){
    const d=panel.data,m=d.metrics,bound=d.bound_points[f].length,load=m.load_nN[f];let values;
    if(stageKey==='g3a'&&panel.role==='mechanism')values=[[m.phase[f],'mechanism state'],[d.time[f].toFixed(1)+' s','time'],[Math.max(...d.protrusion_lengths[f]).toFixed(2)+' µm','protrusion length'],[bound,'bound clutches']];
    else if(stageKey==='g3a')values=[[load.toFixed(3)+' nN','clutch load'],[m.angle_mdeg[f].toFixed(3)+' mdeg','true max Δ angle'],[m.max_bead_dr_um[f].toFixed(5)+' µm','max bead Δr'],[new Set(d.bound_fibre_ids[f]).size,'loaded fibres']];
    else if(stageKey==='g3b'){const angles=d.active_sectors[f].map(s=>(d.cell_angle[f]+2*Math.PI*s/d.config.n_sectors)*180/Math.PI%360);values=[[angles.map(x=>x.toFixed(0)+'°').join(', '),'active θ'],[m.switches[f],'sector switches'],[bound,'bound'],[load.toFixed(3)+' nN','|clutch load|']];}
    else values=[[m.cell_dr_um[f].toFixed(4)+' µm','cell Δr'],[m.cell_rotation_deg[f].toFixed(4)+'°','cell rotation'],[bound,'bound'],[load.toFixed(3)+' nN','|clutch load|']];
    document.querySelector('#'+panel.metricId).innerHTML=values.map(([v,k])=>`<div class="metric"><b>${v}</b><small>${k}</small></div>`).join('');
  }
  function drawPanel(panel,f){
    const d=panel.data,posCount=d.initial_positions.length,attached=new Set(d.bound_fibre_ids[f]),contacts=new Set(d.contact_fibers);
    panel.edgeEls.forEach((node,i)=>{const a=xy(d,f,d.edges[i][0]),b=xy(d,f,d.edges[i][1]),fiber=d.edge_fiber[i];node.setAttribute('x1',a[0]);node.setAttribute('y1',a[1]);node.setAttribute('x2',b[0]);node.setAttribute('y2',b[1]);node.setAttribute('stroke',attached.has(fiber)?'#d45f32':contacts.has(fiber)?'#247c75':'#48636f');node.setAttribute('stroke-width',attached.has(fiber)?.75:contacts.has(fiber)?.48:.33)});
    for(let i=0;i<posCount;i++){const node=panel.beadEls[i],p=xy(d,f,i),fiber=d.bead_fiber[i],color=attached.has(fiber)?'#d45f32':contacts.has(fiber)?'#247c75':'#48636f';node.setAttribute('cx',p[0]);node.setAttribute('cy',p[1]);node.setAttribute('stroke',color)}
    panel.anchors.querySelectorAll('rect').forEach(node=>{const p=xy(d,f,+node.dataset.i);node.setAttribute('x',p[0]-.375);node.setAttribute('y',p[1]-.375)});
    panel.links.replaceChildren();d.crosslinks.forEach(link=>{const a=material(d,f,link,'a'),b=material(d,f,link,'b'),mid=[(a[0]+b[0])/2,(a[1]+b[1])/2];panel.links.append(el('line',{x1:a[0],y1:a[1],x2:b[0],y2:b[1],class:'xconnector'}),el('rect',{x:mid[0]-.34,y:mid[1]-.34,width:.68,height:.68,transform:`rotate(45 ${mid[0]} ${mid[1]})`,class:'xmark'}))});
    const c=d.cell_center[f];panel.cell.setAttribute('cx',c[0]);panel.cell.setAttribute('cy',c[1]);panel.initialCell.style.display=stageKey==='g3c'?'':'none';panel.trail.setAttribute('points',d.cell_center.slice(0,f+1).map(x=>x.join(',')).join(' '));
    panel.activity.replaceChildren();d.activity[f].forEach((a,s)=>{const angle=d.cell_angle[f]+2*Math.PI*s/d.config.n_sectors,r=1.055*d.config.cell_radius;panel.activity.append(el('circle',{cx:c[0]+r*Math.cos(angle),cy:c[1]+r*Math.sin(angle),r:.12+1.05*a/Math.max(...d.activity[f],1e-9),class:'activity'}))});
    panel.protrusions.replaceChildren();const active=new Set(d.active_sectors[f]),attachedSectors=new Set(d.bound_sector_ids[f]);d.protrusion_lengths[f].forEach((length,s)=>{if(length<=.02)return;const angle=d.cell_angle[f]+2*Math.PI*s/d.config.n_sectors,base=[c[0]+d.config.cell_radius*Math.cos(angle),c[1]+d.config.cell_radius*Math.sin(angle)],tip=d.protrusion_tips[f][s],line=el('line',{x1:base[0],y1:base[1],x2:tip[0],y2:tip[1],class:active.has(s)?'protrusion':'retracting'}),dot=el('circle',{cx:tip[0],cy:tip[1],r:active.has(s)?.42:.25,class:'tip'+(attachedSectors.has(s)?' attached':'')});panel.protrusions.append(line,dot)});
    panel.clutches.replaceChildren();d.bound_points[f].forEach((point,i)=>{const motor=d.motor_points[f][i],force=d.clutch_forces_nN[f][i];panel.clutches.append(el('line',{x1:point[0],y1:point[1],x2:motor[0],y2:motor[1],class:'clutch'}),el('circle',{cx:point[0],cy:point[1],r:.36,class:'bound'}),el('line',{x1:point[0],y1:point[1],x2:point[0]+1.8*force[0],y2:point[1]+1.8*force[1],class:'force','marker-end':`url(#arrow-${panel.scene.id})`}))});
    panel.vectors.replaceChildren();const vectorOn=document.querySelector('#vectors').checked&&(panel.role==='response'||stageKey!=='g3a'),scale=+document.querySelector('#vectorScale').value;if(vectorOn){for(let i=0;i<posCount;i+=8){if(d.fixed[i])continue;const a=d.initial_positions[i],b=xy(d,f,i);panel.vectors.append(el('line',{x1:a[0],y1:a[1],x2:a[0]+scale*(b[0]-a[0]),y2:a[1]+scale*(b[1]-a[1]),class:'vector'}))}}
    setView(panel,f);metrics(panel,f);
  }
  function mappedFrame(panel){return Math.round(frame/(frameCount-1)*(panel.data.time.length-1))}
  function draw(){panels.forEach(panel=>drawPanel(panel,mappedFrame(panel)));const reference=panels[0].data,index=mappedFrame(panels[0]);document.querySelector('#time').textContent=reference.time[index].toFixed(1)+' s';document.querySelector('#scrub').value=100*frame/(frameCount-1)}

  document.title=spec.title;document.querySelector('#title').textContent=spec.title;document.querySelector('#subtitle').textContent=spec.subtitle;document.querySelector('#what').textContent=spec.what;document.querySelector('#equation').textContent=spec.equation;document.querySelector('#view').value=spec.view;document.querySelector('#vectors').checked=spec.vectors;document.querySelector('#gifLink').href=spec.gif;
  panels=[setup('left','leftMetrics',spec.left,'#d45f32'),setup('right','rightMetrics',spec.right,'#247c75')];
  const frameCount=Math.max(...panels.map(p=>p.data.time.length));document.querySelector('#scrub').max=100;
  document.querySelector('#play').onclick=()=>{playing=!playing;document.querySelector('#play').textContent=playing?'Pause':'Play'};
  document.querySelector('#step').onclick=()=>{playing=false;frame=Math.min(frame+1,frameCount-1);draw()};
  document.querySelector('#scrub').oninput=e=>{playing=false;frame=Math.round(+e.target.value/100*(frameCount-1));draw()};
  document.querySelector('#view').onchange=draw;document.querySelector('#vectors').onchange=draw;document.querySelector('#vectorScale').onchange=draw;
  function tick(time){if(playing&&time-last>110){frame=(frame+1)%frameCount;last=time;draw()}requestAnimationFrame(tick)}
  draw();requestAnimationFrame(tick);
})();
