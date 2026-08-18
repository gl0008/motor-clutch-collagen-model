(function(){
  const stages={
    g3a:{
      title:'G3A · spatial material-point clutch',
      animation:'figures/g3/g3a_material_point_clutches.gif',
      summary:'figures/g3/g3a_summary.png',
      label:'Fixed-cell single-fibre attachment test',
      watch:'A collagen fibre is a connected bead chain. Red spokes end at bound material points; the attachment follows the same segment coordinate as beads deform.',
      cell:'Fixed center and fixed angle. The single active sector is prescribed only to isolate clutch mechanics.',
      evidence:'Material identity, local force projection, force/torque conservation and rotation-covariance gates pass.',
      metrics:'15 s · FOI 0.6046 → 0.6239 · maximum bead displacement 0.395 µm · maximum 173/200 bound clutches.',
      boundary:'This does not test emergent direction selection or migration. Extending the pull eventually reaches the no-overlap guard.'
    },
    g3b:{
      title:'G3B · collagen-guided protrusions',
      animation:'figures/g3/g3b_emergent_protrusions.gif',
      summary:'figures/g3/g3b_summary.png',
      label:'Fixed-cell aligned-collagen guidance test',
      watch:'The two red active-sector rings turn over among 24 initially unbiased probes. Larger/brighter sector points mean more nearby collagen aligned with that outward direction.',
      cell:'Fixed center and fixed angle so protrusion direction is not confused with cell translation.',
      evidence:'Aligned guidance, ± director symmetry, 30° rotation covariance, feedback ablation and no-fibre controls pass. The overall calibration gate remains blocked.',
      metrics:'Aligned 20-seed mean nematic guidance 0.726 (95% CI 0.626–0.819) · 55/45 axis-sign split · 3.03° rotation error.',
      boundary:'Seven of twenty random-isotropic runs pulled collagen inside the rigid cell. Until an explicit contact law is added and all controls are rerun, G3B is not formally validated.'
    },
    g3c:{
      title:'G3C · reaction-driven translation and rotation',
      animation:'figures/g3/g3c_translation_rotation.gif',
      secondary:'figures/g3/g3c_direction_controls.gif',
      summary:'figures/g3/g3c_summary.png',
      label:'Released rigid-cell force/torque smoke test',
      watch:'Spatial clutch vectors create an equal-and-opposite blue reaction arrow and torque. The blue trajectory is the resulting cell path, not a prescribed velocity.',
      cell:'Both x/y translation and body rotation are released under linear overdamped drag.',
      evidence:'The reaction mechanism and hidden-drive audit exist, but the formal G3C campaign is halted behind the failed upstream G3B geometry gate.',
      metrics:'30 s asymmetric fixture · net displacement 0.00187 µm · final rotation 6.68×10⁻⁵ rad. Values are not display-amplified.',
      boundary:'This tiny motion is a mechanism demonstration, not realistic tumor migration. The partial 173-run calibration data are retained for audit and excluded from evidence.'
    }
  };
  const key=new URLSearchParams(location.search).get('stage');
  const stage=stages[key] ? key : 'g3a';
  const data=stages[stage];
  const byId=id=>document.getElementById(id);
  document.title=data.title;
  byId('title').textContent=data.title;
  byId('watch').textContent=data.watch;
  byId('cellState').textContent=data.cell;
  byId('evidence').textContent=data.evidence;
  byId('animationLabel').textContent=data.label;
  byId('metrics').textContent=data.metrics;
  byId('boundary').textContent=data.boundary;
  byId('animation').src=data.animation;
  byId('summaryImage').src=data.summary;
  document.querySelectorAll('nav a').forEach(link=>link.classList.toggle('current',link.dataset.stage===stage));
  const restart=(image,path)=>{image.src='';requestAnimationFrame(()=>{image.src=path+'?restart='+Date.now();});};
  byId('restart').addEventListener('click',()=>restart(byId('animation'),data.animation));
  if(data.secondary){
    byId('secondary').hidden=false;
    byId('secondaryAnimation').src=data.secondary;
    byId('restartSecondary').addEventListener('click',()=>restart(byId('secondaryAnimation'),data.secondary));
  }
})();
