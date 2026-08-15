const T=String.raw;
const REPO='https://github.com/gl0008/motor-clutch-collagen-model/blob/main/';
const PAPERS={
  lee:{title:'Lee et al. (2014) · 3D computational collagen network mechanics',href:'https://doi.org/10.1371/journal.pone.0111896'},
  abhilash:{title:'Abhilash et al. (2014) · remodeling of fibrous ECM by contractile cells',href:'https://doi.org/10.1016/j.bpj.2014.08.029'},
  wang:{title:'Wang et al. (2014) · tension-driven long-range force transmission',href:'https://doi.org/10.1016/j.bpj.2014.09.044'},
  notbohm:{title:'Notbohm et al. (2015) · microbuckling in fibrous matrices',href:'https://arxiv.org/abs/1407.3510'},
  bell:{title:'Bell (1978) · force-dependent adhesion kinetics',href:'https://doi.org/10.1126/science.347575'},
  bangasser:{title:'Bangasser & Odde (2013) · motor–clutch master equation',href:'https://doi.org/10.1007/s12195-013-0296-5'},
  prahl:{title:'Prahl et al. (2020) · opposing motor–clutch modules in 1D migration',href:'https://doi.org/10.1016/j.bpj.2020.01.048'},
  steinwachs:{title:'Steinwachs et al. (2016) · forces in collagen networks',href:'https://doi.org/10.1038/nmeth.3685'},
  speed:{title:'Aguilar-Rojas et al. (2024) · MDA-MB-231 speeds in collagen',href:'https://pmc.ncbi.nlm.nih.gov/articles/PMC11327030/'},
  han:{title:'Han et al. (2018) · long-ranged matrix effects from cell contraction',href:'https://pmc.ncbi.nlm.nih.gov/articles/PMC5910866/'},
  ban:{title:'Ban et al. (2018) · plastic deformation in collagen networks',href:'https://doi.org/10.1016/j.bpj.2017.11.3739'},
  kim:{title:'Kim et al. (2017) · stress-induced plasticity of dynamic collagen',href:'https://doi.org/10.1038/s41467-017-01011-7'},
  adebowale:{title:'Adebowale et al. (2021) · stress relaxation and migration',href:'https://doi.org/10.1038/s41563-021-00981-w'}
};
const paper=(key,support,boundary)=>({...PAPERS[key],support,boundary});

const G1_V2_EQUATIONS=[
  {title:'Bead–spring fibre mechanics',intro:'Each visible fibre is a connected bead chain. The initially curved geometry is stress-free.',equations:[
    {name:'Bond length and extension',tex:T`\[\ell_{ij}=\|\mathbf r_j-\mathbf r_i\|,\qquad \Delta\ell_{ij}=\ell_{ij}-\ell^0_{ij}\]`,explain:'The bond only reacts to change from its own initial length.',units:'µm for length; G1 force scale is arbitrary.'},
    {name:'Stretching energy',tex:T`\[U_s=\frac12 k_s\sum_{(i,j)}(\Delta\ell_{ij})^2\]`,explain:'Hookean axial energy stored by every bead-to-bead spring.'},
    {name:'Stretching force',tex:T`\[\mathbf F^{s}_{i\leftarrow j}=k_s\Delta\ell_{ij}\frac{\mathbf r_j-\mathbf r_i}{\ell_{ij}},\qquad \mathbf F^{s}_{j\leftarrow i}=-\mathbf F^{s}_{i\leftarrow j}\]`,explain:'Equal and opposite forces are applied to the two beads.'},
    {name:'Reference-curvature change',tex:T`\[\mathbf c_i=(\mathbf r_{i-1}-2\mathbf r_i+\mathbf r_{i+1})-\mathbf c_i^0\]`,explain:'Subtracting the initial curvature prevents the generated wavy fibre from straightening before loading.'},
    {name:'Bending energy and bead forces',tex:T`\[U_b=\frac12k_b\sum_i\|\mathbf c_i\|^2,\quad(\mathbf F^b_{i-1},\mathbf F^b_i,\mathbf F^b_{i+1})=(-k_b\mathbf c_i,2k_b\mathbf c_i,-k_b\mathbf c_i)\]`,explain:'A discrete curvature penalty makes a bead chain act like a flexible fibre rather than disconnected axial springs.'}
  ]},
  {title:'Permanent material-point crosslinks',equations:[
    {name:'Point on a bond',tex:T`\[\mathbf p(\alpha)=(1-\alpha)\mathbf r_i+\alpha\mathbf r_j,\qquad 0\le\alpha\le1\]`,explain:'A crosslink is attached to a material position along a segment, not necessarily at a bead.'},
    {name:'Crosslink deformation',tex:T`\[\boldsymbol\delta_x=(\mathbf p_b-\mathbf p_a)-\boldsymbol\ell_x^0\]`,explain:'The link is stress-free at its initial rest vector.'},
    {name:'Crosslink energy and force',tex:T`\[U_x=\frac12k_x\|\boldsymbol\delta_x\|^2,\qquad \mathbf F_x=k_x\boldsymbol\delta_x\]`,explain:'This transfers force between different fibres while allowing free rotation at the crossing.',key:true},
    {name:'Force distribution to segment beads',tex:T`\[(\mathbf F_i,\mathbf F_j)=((1-\alpha)\mathbf F_x,\alpha\mathbf F_x)\]`,explain:'Linear shape-function weights preserve the crosslink force and moment along the segment.'}
  ]},
  {title:'Hybrid contact + Gaussian load distribution',equations:[
    {name:'Closest point on a segment',tex:T`\[\alpha^*=\operatorname{clip}\!\left(\frac{(\mathbf r_c-\mathbf a)\cdot(\mathbf b-\mathbf a)}{\|\mathbf b-\mathbf a\|^2},0,1\right),\qquad \mathbf p=\mathbf a+\alpha^*(\mathbf b-\mathbf a)\]`,explain:'Only the closest material point on each fibre is considered.'},
    {name:'Hard contact eligibility',tex:T`\[0\le d_s=\|\mathbf p-\mathbf r_c\|-R\le d_{contact},\qquad |\Delta\theta|\le\theta_{half}\]`,explain:'The Gaussian is not applied to the entire network. A point must first lie in the protrusion contact shell.',key:true},
    {name:'Normalized Gaussian weights',tex:T`\[\tilde w_m=e^{-d_{s,m}^2/\sigma_c^2},\qquad w_m=\frac{\tilde w_m}{\sum_n\tilde w_n},\qquad \sum_mw_m=1\]`,explain:'Closer eligible fibres receive more of the prescribed total force, without limiting coupling to discrete nearest beads.',key:true},
    {name:'Patch and nodal forces',tex:T`\[\mathbf f_m=F_{tot}w_m\mathbf n_{in,m},\qquad(\mathbf F_i,\mathbf F_j)=((1-\alpha_m)\mathbf f_m,\alpha_m\mathbf f_m)\]`,explain:'The total pull is conserved and spread continuously along the contacted segment.'},
    {name:'Rigid-cell repulsion',tex:T`\[\mathbf F_i^{rep}=k_{rep}\max(0,R+\epsilon-\|\mathbf r_i-\mathbf r_c\|)\frac{\mathbf r_i-\mathbf r_c}{\|\mathbf r_i-\mathbf r_c\|}\]`,explain:'A soft penalty prevents collagen beads from passing through the rigid cell.'}
  ]},
  {title:'Overdamped motion, loading and observables',equations:[
    {name:'Overdamped bead force balance',tex:T`\[\zeta_b\dot{\mathbf r}_i=\mathbf F_i^s+\mathbf F_i^b+\mathbf F_i^x+\mathbf F_i^{rep}+\mathbf F_i^{act}\]`,explain:'Velocity is proportional to the instantaneous sum of elastic and active forces; inertia is omitted.',key:true},
    {name:'Explicit position update',tex:T`\[\mathbf r_i^{n+1}=\mathbf r_i^n+\Delta t\,\frac{\mathbf F_i^{tot,n}}{\zeta_b}\]`,explain:'Free beads take an Euler step; every fibre end in G1 is then reset to its initial position.'},
    {name:'Force ramp',tex:T`\[F_{tot}(t)=F_{max}\min(1,t/t_{ramp})\]`,explain:'The pull ramps on to reduce a numerical impulse.'},
    {name:'Radial alignment order',tex:T`\[S_r=\left\langle2(\hat{\mathbf t}\cdot\hat{\mathbf e}_r)^2-1\right\rangle\]`,explain:'Sᵣ=1 is radial, 0 is isotropic in 2D, and −1 is tangential.',key:true},
    {name:'RMS network displacement',tex:T`\[u_{RMS}=\sqrt{\frac1{N_{free}}\sum_{i\in free}\|\mathbf r_i-\mathbf r_i^0\|^2}\]`,explain:'A scalar summary of how far mobile beads move.'},
    {name:'Power input and drag dissipation',tex:T`\[P_{act}=\sum_i\mathbf F_i^{act}\cdot\dot{\mathbf r}_i,\qquad P_{drag}=\zeta_b\sum_i\|\dot{\mathbf r}_i\|^2\]`,explain:'These distinguish energy supplied by the cell from viscous dissipation.'}
  ]}
];

const G2_BASE_EQUATIONS=[
  {title:'Physical fibre calibration',equations:[
    {name:'Cross-section and second moment',tex:T`\[A=\frac{\pi d^2}{4},\qquad I=\frac{\pi d^4}{64}\]`,explain:'A circular effective fibre diameter converts modulus into axial and bending rigidities.',units:'d=0.30 µm in the baseline.'},
    {name:'Axial and bending rigidity',tex:T`\[EA=E A,\qquad EI=E I\]`,explain:'The baseline uses E=32 MPa and converts 1 MPa to 1000 nN/µm².'},
    {name:'Resolution-aware coefficients',tex:T`\[k_t=\frac{EA}{\ell_0},\qquad k_c=\rho_c k_t,\qquad k_b=\frac{EI}{\ell_0^3}\]`,explain:'Tension uses kₜ; compression uses ρc=0.10 of it to represent microbuckling; bending scales with bead spacing.',key:true}
  ]},
  {title:'Tension/compression and bending',equations:[
    {name:'Bond strain and branch stiffness',tex:T`\[\varepsilon_{ij}=\frac{\ell_{ij}-\ell^0_{ij}}{\ell^0_{ij}},\qquad k(\Delta\ell)=\begin{cases}k_t,&\Delta\ell\ge0\\k_c,&\Delta\ell<0\end{cases}\]`,explain:'The sign of extension selects the tensile or compression-softened branch.'},
    {name:'Axial force and energy',tex:T`\[\mathbf F^{s}_{i\leftarrow j}=k(\Delta\ell)\Delta\ell\frac{\mathbf r_j-\mathbf r_i}{\ell_{ij}},\qquad U_s=\frac12\sum k(\Delta\ell)(\Delta\ell)^2\]`,explain:'Connected beads visibly trace the fibre while the spring stores axial energy.'},
    {name:'Reference-curvature mechanics',tex:T`\[\mathbf c_i=(\mathbf r_{i-1}-2\mathbf r_i+\mathbf r_{i+1})-\mathbf c_i^0,\quad U_b=\frac12k_b\sum_i\|\mathbf c_i\|^2\]`,explain:'The generated curved fibre is stress-free; only changes in curvature cost energy.'},
    {name:'Bending bead forces',tex:T`\[(\mathbf F^b_{i-1},\mathbf F^b_i,\mathbf F^b_{i+1})=(-k_b\mathbf c_i,2k_b\mathbf c_i,-k_b\mathbf c_i)\]`,explain:'This is the negative gradient of the discrete curvature energy.'}
  ]},
  {title:'Freely hinged inter-fibre crosslinks',equations:[
    {name:'Material points',tex:T`\[\mathbf p_a=(1-\alpha_a)\mathbf r_i+\alpha_a\mathbf r_j,\qquad \mathbf p_b=(1-\alpha_b)\mathbf r_k+\alpha_b\mathbf r_l\]`,explain:'The two connected points may lie between beads.'},
    {name:'Link deformation, force and energy',tex:T`\[\boldsymbol\delta_x=(\mathbf p_b-\mathbf p_a)-\boldsymbol\ell_x^0,\quad \mathbf F_x=k_x\boldsymbol\delta_x,\quad U_x=\frac12k_x\|\boldsymbol\delta_x\|^2\]`,explain:'The link resists separation but imposes no preferred crossing angle.',key:true},
    {name:'Nodal force transfer',tex:T`\[(\mathbf F_i,\mathbf F_j,\mathbf F_k,\mathbf F_l)=((1-\alpha_a)\mathbf F_x,\alpha_a\mathbf F_x,-(1-\alpha_b)\mathbf F_x,-\alpha_b\mathbf F_x)\]`,explain:'This is the mechanism that allows a directly pulled fibre to move indirectly connected fibres.'},
    {name:'Boundary-percolation acceptance gate',tex:T`\[f_{connected}=\frac{N_{fibres\ connected\ to\ boundary}}{N_{fibres}}\ge0.85,\qquad \mathcal C_{contact}\subseteq\mathcal C_{boundary}\]`,explain:'A generated network is rejected unless all contact fibres and at least 85% of fibres connect to the anchored outer boundary.',key:true}
  ]},
  {title:'Contact-region + Gaussian coupling',equations:[
    {name:'Closest point on each fibre',tex:T`\[\alpha^*=\operatorname{clip}\!\left(\frac{(\mathbf r_c-\mathbf a)\cdot(\mathbf b-\mathbf a)}{\|\mathbf b-\mathbf a\|^2},0,1\right),\quad \mathbf p=\mathbf a+\alpha^*(\mathbf b-\mathbf a)\]`,explain:'The coupling follows the continuous bead–spring segment rather than snapping only to bead centers.'},
    {name:'Hard contact shell and protrusion sector',tex:T`\[0\le\|\mathbf p-\mathbf r_c\|-R\le3\ \mu m,\qquad |\Delta\theta|\le30^\circ\]`,explain:'Only physical near-surface candidates may receive direct cell force.',key:true},
    {name:'Gaussian weights among eligible contacts',tex:T`\[w_m=\frac{\exp(-d_{s,m}^2/\sigma_c^2)}{\sum_n\exp(-d_{s,n}^2/\sigma_c^2)},\qquad \sigma_c=1.5\ \mu m\]`,explain:'Closer contacts receive more force, but ∑w=1 keeps the total force fixed.',key:true},
    {name:'Patch and bead forces',tex:T`\[\mathbf f_m=F_{tot}w_m\mathbf n_{in,m},\qquad(\mathbf F_i,\mathbf F_j)=((1-\alpha_m)\mathbf f_m,\alpha_m\mathbf f_m)\]`,explain:'The arrow direction is normal to the cell surface; the display arrow is now shortened without changing this physics.'},
    {name:'No-penetration penalty',tex:T`\[\mathbf F_i^{rep}=k_{rep}\max(0,R+\epsilon-r_i)\hat{\mathbf e}_{r,i}\]`,explain:'This maintains a rigid-cell exclusion region.'}
  ]},
  {title:'Overdamped dynamics, boundary condition and validation observables',equations:[
    {name:'Overdamped bead force balance',tex:T`\[\boxed{\zeta_b\dot{\mathbf r}_i=\mathbf F_i^s+\mathbf F_i^b+\mathbf F_i^x+\mathbf F_i^{rep}+\mathbf F_i^{act}}\]`,explain:'This is the central equation: stored elastic forces plus active pull determine velocity against drag. There is no bead inertia and no SLS state.',key:true},
    {name:'Time integration',tex:T`\[\mathbf r_i^{n+1}=\mathbf r_i^n+\Delta t\,\mathbf F_i^{tot,n}/\zeta_b\]`,explain:'Only non-boundary beads update; beads inside the 2.5 µm outer boundary band remain fixed.'},
    {name:'Force ramp',tex:T`\[F_{tot}(t)=F_{max}\min(1,t/t_{ramp}),\qquad F_{max}=5\ \mathrm{nN}\]`,explain:'The accepted elastic baseline ramps to a 5 nN total pull.'},
    {name:'Radial alignment by distance shell',tex:T`\[S_r^{(q)}=\left\langle2(\hat{\mathbf t}\cdot\hat{\mathbf e}_r)^2-1\right\rangle_{q},\quad q\in\{0\!:\!10,10\!:\!25,25\!:\!50\}\ \mu m\]`,explain:'Near, intermediate and far shells distinguish local alignment from crosslink-mediated transmission.',key:true},
    {name:'Shell displacement and RMS displacement',tex:T`\[\bar u_q=\langle\|\mathbf r_i-\mathbf r_i^0\|\rangle_q,\qquad u_{RMS}=\sqrt{N_{free}^{-1}\sum_{free}\|\mathbf r_i-\mathbf r_i^0\|^2}\]`,explain:'These quantify motion rather than relying on a visually dramatic animation.'},
    {name:'Input and dissipation',tex:T`\[P_{act}=\sum_i\mathbf F_i^{act}\cdot\dot{\mathbf r}_i,\qquad P_{drag}=\zeta_b\sum_i\|\dot{\mathbf r}_i\|^2\]`,explain:'Active work is dissipated by drag or stored transiently in stretch, bend and links.'},
    {name:'Energy-balance check',tex:T`\[\Delta(U_s+U_b+U_x)\approx\int(P_{act}-P_{drag})\,dt\]`,explain:'A small relative mismatch checks that the force and integration bookkeeping are consistent.'},
    {name:'Minimum cell gap',tex:T`\[g_{min}=\min_i\|\mathbf r_i-\mathbf r_c\|-R\ge0\]`,explain:'A nonnegative gap verifies that fibres do not pass through the cell.'}
  ]}
];

window.NOTEBOOK_MODELS={};

NOTEBOOK_MODELS['g1-v0']={
  generation:'G1',version:'V0',status:'question-forming prototype',title:'SLS single protrusion and triangular lattice',
  subtitle:'The project began by asking whether matrix stress relaxation changes clutch loading and cluster lifetime before trying to simulate a collagen network.',
  question:'If two substrates have the same stiffnesses but different SLS relaxation times, does faster relaxation reduce clutch loading rate and extend adhesion-cluster lifetime?',
  units:'The 1D model is nondimensional. The triangular lattice uses prototype force/length/time scales rather than calibrated nN–µm collagen parameters.',
  purpose:'把因果鏈先縮成一隻 protrusion：actin motion → clutch loading → SLS relaxation → force-dependent rupture → cluster lifetime；之後才嘗試 2D bead network。',
  brief:[
    '「先建立固定 cell body、單一 protrusion 的 1D 模型；第一階段 endpoint 是 clutch-cluster lifetime。」',
    '「把原本 bead–spring bond 的 spring 換成 SLS，fast/slow relaxation 使用相同 stochastic stream。」',
    '「10 beads、9 bonds；每條 bond stiffness 乘 9，讓整條 chain stiffness 不隨 bead resolution 改變。」'
  ],
  chain:['actin moves','clutches load','SLS relaxes','Bell rupture','cluster lifetime'],
  artifacts:[
    {title:'1D single-protrusion source',href:REPO+'versions/v0_sls_prototypes/single_protrusion_1d/single_protrusion.py',note:'Fast/slow SLS, 12 stochastic clutches, motor force–velocity law and cluster episodes.'},
    {title:'2D triangular-lattice source',href:REPO+'versions/v0_sls_prototypes/triangular_lattice_2d/network.py',note:'Diluted triangular network with SLS axial segments, elastic bending and overdamped bead motion.'}
  ],
  source:REPO+'versions/v0_sls_prototypes/README.md',
  equationGroups:[
    {title:'Nondimensional scales and serial-chain resolution',equations:[
      {name:'Reference scales',tex:T`\[F_*=F_b,\qquad K_*=k_c,\qquad L_*=\frac{F_b}{k_c},\qquad t_{load}=\frac{F_b}{k_cv_0}\]`,explain:'Force is measured by the Bell rupture scale, stiffness by clutch stiffness, extension by force/stiffness and time by the unloaded time needed to load one rupture force.',key:true},
      {name:'Nondimensional variables',tex:T`\[\tilde F=F/F_*,\quad \tilde x=x/L_*,\quad \tilde t=t/t_{load},\quad \tilde k=k/K_*\]`,explain:'This makes the first experiment mechanism-first rather than tied to uncertain physical calibration.'},
      {name:'Springs in series',tex:T`\[K_{chain}^{-1}=\sum_{b=1}^{N_b}k_{bond}^{-1}=\frac{N_b}{k_{bond}}\quad\Rightarrow\quad k_{bond}=N_bK_{chain}\]`,explain:'With 9 identical bonds, each bond must be 9 times the desired chain stiffness.',key:true},
      {name:'Both SLS stiffness limits scale the same way',tex:T`\[k_{0,bond}=N_bK_{0,chain},\qquad k_{\infty,bond}=N_bK_{\infty,chain}\]`,explain:'Therefore the instantaneous and long-time end-to-end stiffnesses remain resolution independent.'}
    ]},
    {title:'Standard linear solid constitutive law',equations:[
      {name:'Parallel equilibrium spring + Maxwell force',tex:T`\[F=K_\infty x+q,\qquad K_m=K_0-K_\infty\]`,explain:'K∞ stores the long-time force; q is the relaxing Maxwell-branch force.',key:true},
      {name:'Maxwell-state evolution',tex:T`\[\dot q=K_m\dot x-\frac{q}{\tau},\qquad \eta=\tau K_m\]`,explain:'Extension rate loads the Maxwell spring while the dashpot relaxes q.',key:true},
      {name:'Exact finite-step update',tex:T`\[q_{n+1}=e^{-\Delta t/\tau}q_n+K_m\tau(1-e^{-\Delta t/\tau})\frac{x_{n+1}-x_n}{\Delta t}\]`,explain:'The code integrates q exactly when extension changes linearly over one time step.'},
      {name:'Step-strain relaxation check',tex:T`\[F(t)=K_\infty\Delta x+(K_0-K_\infty)\Delta x\,e^{-t/\tau}\]`,explain:'This analytic curve is the unit test for a single SLS element.'}
    ]},
    {title:'1D motor–clutch coupling',equations:[
      {name:'Tension-only clutch spring',tex:T`\[F_i=k_c\max(x_a-x_{bind,i}-x_s,0)\]`,explain:'A bound clutch carries tension only; xₛ is the substrate-chain extension.'},
      {name:'Quasi-static endpoint balance',tex:T`\[F_{SLS}(x_s,q)=\sum_{i\in bound}F_i\]`,explain:'Massless internal beads make every serial bond carry the same traction.',key:true},
      {name:'Motor force–velocity relation',tex:T`\[v_a=v_0\max\!\left(0,1-\frac{\sum_iF_i}{F_{stall}}\right),\qquad x_a^{n+1}=x_a^n+v_a\Delta t\]`,explain:'Actin slows as total clutch traction approaches stall force.',key:true},
      {name:'Bell slip-bond off-rate',tex:T`\[k_{off,i}=k_{off}^0\exp(F_i/F_b)\]`,explain:'Higher clutch force accelerates rupture.'},
      {name:'Discrete event probabilities',tex:T`\[P_{on}=1-e^{-k_{on}\Delta t},\qquad P_{off,i}=1-e^{-k_{off,i}\Delta t}\]`,explain:'Counter-addressed uniform random numbers make fast and slow cases share proposed events.'},
      {name:'Cluster lifetime',tex:T`\[T_{cluster}=t(N_{bound}\to0)-t(0\to N_{bound}>0)\]`,explain:'An episode starts at the first binding and ends when the bound count returns to zero.'},
      {name:'Force-loading rate',tex:T`\[\dot F_T\approx\frac{F_T^{n+1}-F_T^n}{\Delta t},\qquad F_T=\sum_iF_i\]`,explain:'The central mechanistic expectation was that faster relaxation lowers positive loading rate.'}
    ]},
    {title:'Triangular-lattice extension',equations:[
      {name:'SLS axial force per edge',tex:T`\[F_e=K_\infty(\ell_e-\ell_e^0)+q_e,\qquad \dot q_e=(K_0-K_\infty)\dot\ell_e-q_e/\tau\]`,explain:'Every triangular-lattice edge receives its own SLS internal state.'},
      {name:'Elastic bending',tex:T`\[\mathbf c_i=\mathbf r_{i-1}-2\mathbf r_i+\mathbf r_{i+1},\qquad U_b=\frac12k_b\sum_i\|\mathbf c_i\|^2\]`,explain:'Bending applies only to collinear triplets within a lattice fibre.'},
      {name:'Overdamped bead motion',tex:T`\[\zeta_b\dot{\mathbf r}_i=\mathbf F_i^{SLS}+\mathbf F_i^b+\mathbf F_i^{clutch}\]`,explain:'This combines material relaxation inside bonds with external drag on bead motion; they are different dissipative mechanisms.'},
      {name:'Point-cell translation',tex:T`\[\zeta_c\dot{\mathbf r}_c=-\sum_i\mathbf F_i^{clutch}\]`,explain:'The reaction to clutch forces moves the point cell.'}
    ]}
  ],
  algorithm:[
    'Create 12 unbound clutches, an unloaded SLS chain and a fixed cell body.',
    'Use the same counter-addressed random number for the same step/clutch/event in fast and slow conditions.',
    'Update Bell rupture and new binding events.',
    'Update actin velocity from the previous total traction.',
    'Solve the tension-only active set and SLS/clutch force balance for the new chain extension.',
    'Record loading, relaxation, rupture, cluster failure and lifetime; repeat for τ=0.1 and τ=10.',
    'The triangular extension replaces the serial chain with a diluted lattice and advances all free beads overdamped.'
  ],
  papers:[
    paper('adebowale','Motivates the hypothesis that stress-relaxing substrates can change protrusion/migration behavior and therefore clutch loading history.','That a collagen bead bond must be an SLS element, or that these dimensionless parameters are calibrated collagen values.'),
    paper('bell','Supports the exponential force-accelerated slip-bond off-rate used for clutch rupture.','That every integrin adhesion behaves as a single Bell slip bond.'),
    paper('bangasser','Supports stochastic effective clutches, actin loading and force–velocity feedback.','The exact 12-clutch coarse graining or the selected nondimensional default values.'),
    paper('lee','Supports representing collagen-like structures with connected beads and springs.','The triangular lattice geometry or SLS constitutive choice in V0.')
  ],
  results:[
    {name:'Primary endpoint',value:'Distribution of cluster lifetimes from first binding until complete cluster failure.'},
    {name:'Paired comparison',value:'Fast and slow relaxation share the same event-addressed random stream until mechanics changes their hazards.'},
    {name:'Resolution test',value:'5, 9 and 19 serial bonds reproduce the same end-to-end K₀, K∞ and τ.'},
    {name:'Interpretation status',value:'Mechanism prototype only; not a calibrated collagen network.'}
  ],
  limitation:{heading:'HOWEVER · why V0 could not be the collagen baseline',text:'SLS relaxation inside a bond and overdamped bead drag are not the same physics, but V0 combined them before the geometric collagen-network response was understood. The triangular lattice also looked unlike finite collagen fibres.',items:['The 1D chain cannot show crosslink-mediated alignment or spatial force decay.','The triangular lattice hard-codes a regular topology and perimeter rather than finite 20–80 µm fibres.','SLS parameters were dimensionless and not calibrated to collagen.','Too many mechanisms changed simultaneously to attribute cause and effect.'],leadsTo:'G1 V1: temporarily remove SLS and build the simplest few-fibre overdamped elastic baseline.'},
  talk:{script:'“V0 was my hypothesis-forming stage. I isolated motor loading, SLS relaxation and Bell rupture in one dimension, then tested whether the idea survived a bead network. It taught me the clutch equations and resolution scaling, but it did not yet represent collagen geometry, so I did not carry SLS into the next baseline.”',assumptions:['Quasi-static axial equilibrium in the 1D chain.','Tension-only clutches and Bell slip-bond rupture.','The same K₀, K∞ and τ describe every serial bond after resolution scaling.'],nextChecks:['Calibrate any future SLS term independently from bead drag.','First validate elastic crosslinked fibre response and spatial alignment.','Only reintroduce stress relaxation after the elastic baseline is identifiable.']}
};

NOTEBOOK_MODELS['g1-v1']={
  generation:'G1',version:'V1',status:'minimal few-fibre baseline',title:'Fixed cell pulling five bead–spring fibres',
  subtitle:'A deliberately small model to see stretching, bending, active pulling and overdamped drag before adding crosslinks.',
  question:'Can a fixed circular cell deform a few connected bead chains under an overdamped force balance, and can every force term be visually traced?',
  units:'Geometry is in µm and time in seconds; force and stiffness remain prototype/dimensionless scales.',
  purpose:'先確定 beads 真的連成 fibres，而且受力後是 stretching + bending + drag 的連續變形，不先加入 crosslinks、breaking 或 cell movement。',
  brief:['「先用 few fibres，固定 ball；拉 closest beads，方向 normal to the cell surface。」','「SLS 先拿掉，先看 overdamped stretching、bending、active pull、drag 是否產生 reorganization。」','「beads 應該看起來是一條連起來的線，不是一堆點 pop up。」'],
  chain:['fixed cell','all-bead Gaussian pull','fibre stretch/bend','overdamped motion','visible deformation'],
  simulation:{src:'versions/v1_few_fiber/demo/index.html',note:'This original animation is preserved exactly as the first clearly bead-based visual baseline.'},
  source:REPO+'versions/v1_few_fiber/model.py',
  equationGroups:[
    {title:'Five connected fibres',equations:[
      {name:'Stretching energy and force',tex:T`\[U_s=\frac12k_s\sum_e(\ell_e-\ell_e^0)^2,\qquad \mathbf F_e=k_s(\ell_e-\ell_e^0)\hat{\mathbf t}_e\]`,explain:'Each bead-to-bead edge is an elastic axial spring.'},
      {name:'Reference-curvature bending',tex:T`\[\mathbf c_i=(\mathbf r_{i-1}-2\mathbf r_i+\mathbf r_{i+1})-\mathbf c_i^0,\qquad U_b=\frac12k_b\sum_i\|\mathbf c_i\|^2\]`,explain:'The fibre resists changes from its initial curve.'},
      {name:'All-bead Gaussian weighting',tex:T`\[\tilde w_i=\exp\!\left[-\left(\frac{\|\mathbf r_i-\mathbf r_c\|-R}{\sigma}\right)^2\right],\qquad w_i=\frac{\tilde w_i}{\sum_j\tilde w_j}\]`,explain:'Every bead can receive nonzero weight, even when it is not truly in contact. This is the key flaw later corrected.',key:true},
      {name:'Radial active pull',tex:T`\[\mathbf F_i^{act}=F_{tot}(t)w_i\frac{\mathbf r_c-\mathbf r_i}{\|\mathbf r_c-\mathbf r_i\|}\]`,explain:'The cell pulls beads inward along the local cell normal.'},
      {name:'Force ramp',tex:T`\[F_{tot}(t)=F_{max}\min(1,t/t_{ramp})\]`,explain:'The applied force rises smoothly from zero.'},
      {name:'Overdamped dynamics',tex:T`\[\boxed{\zeta_b\dot{\mathbf r}_i=\mathbf F_i^s+\mathbf F_i^b+\mathbf F_i^{act}}\]`,explain:'With no inertia, the net force directly determines bead velocity.',key:true},
      {name:'Euler update and fixed ends',tex:T`\[\mathbf r_i^{n+1}=\mathbf r_i^n+\Delta t\,\mathbf F_i^{tot}/\zeta_b,\qquad \mathbf r_{end}(t)=\mathbf r_{end}(0)\]`,explain:'Both ends of each finite fibre are fixed.'}
    ]}
  ],
  algorithm:['Generate five gently curved bead chains around a fixed circular cell.','Store each bond length and each triplet curvature as the stress-free reference.','Compute Gaussian weights over all beads and ramp the total inward pulling force.','Sum stretching, bending and active forces.','Divide by bead drag, update free beads and reset fibre ends to their initial locations.','Animate connected springs and beads at every saved time.'],
  papers:[paper('lee','Supports a connected bead–spring representation of collagen fibres.','The five-fibre topology, all-bead Gaussian load or prototype stiffnesses.'),paper('abhilash','Supports examining fibre reorientation under local cellular contraction.','That five isolated fibres are sufficient for long-range remodeling.'),paper('steinwachs','Motivates a breast-cancer-cell/collagen force context.','The uncalibrated force unit used in G1 V1.')],
  results:[{name:'Visual success',value:'Beads remain connected as continuous fibres and deform smoothly rather than appearing as disconnected points.'},{name:'Mechanics included',value:'Stretching, bending, active pull and overdamped drag.'},{name:'Mechanics excluded',value:'Crosslinks, clutch stochasticity, cell translation and permanent remodeling.'},{name:'Main diagnostic',value:'Whether the pulled fibre visibly bends/reorients without numerical instability.'}],
  limitation:{heading:'HOWEVER · the Gaussian was applied too broadly',text:'The all-bead Gaussian gives mathematically nonzero force to remote beads and the five fibres cannot transmit force through crosslinks.',items:['No hard contact region before Gaussian weighting.','No inter-fibre crosslinks, so neighboring fibres cannot move indirectly.','Both ends of every fibre are fixed.','Force/stiffness are not calibrated physical units.'],leadsTo:'G1 V2: add a hard contact shell, material-point Gaussian weights and permanent crosslinks.'},
  talk:{script:'“V1 is my cleanest mechanics sanity check. A bead chain stores stretch and bend energy; the summed force divided by drag gives bead velocity. I fixed the cell and fibre ends so I could see one cause at a time. Its failure is equally important: remote beads still received tiny Gaussian loads and there was no crosslink pathway.”',assumptions:['Five fibres are enough to debug mechanics, not to represent bulk ECM.','Initial curvature is stress-free.','Every fibre end is an anchor.'],nextChecks:['Add contact eligibility before Gaussian weighting.','Connect different fibres and compare with/without crosslink force.','Measure alignment and displacement rather than judging only by eye.']}
};

NOTEBOOK_MODELS['g1-v2']={
  generation:'G1',version:'V2',status:'first crosslinked elastic network',title:'Fixed cell pulling a permanently crosslinked network',
  subtitle:'The first hybrid contact-region + Gaussian model, created to test whether permanent crosslinks spread a local pull into neighboring fibres.',
  question:'Compared with the same fibres without crosslink force, do permanent intersection links transmit displacement and alignment away from direct cell contacts?',
  units:'Geometry is in µm and time in seconds; force remains an arbitrary prototype unit despite a 5-unit total pull.',
  purpose:'fix cell-versus-collagen scale, add permanent crosslinks, then verify overdamped stretching + bending + pulling + drag before allowing any link to break。',
  brief:['「nearest/contact region + Gaussian weighting 的 hybrid。」','「加 crosslinks 來代表力傳到其他 fibres；到這裡 crosslink 都還是不會斷。」','「near fibres fully realign, linked fibres partially remodel, far fibres should barely change; ends fixed for far-field anchoring。」'],
  chain:['hard contact shell','Gaussian share of fixed total force','direct fibres deform','crosslinks transmit force','near/intermediate/far alignment'],
  simulation:{src:'versions/v2_crosslinked_elastic/demo/index.html',note:'The two panels compare the same network without and with permanent crosslink forces.'},
  source:REPO+'versions/v2_crosslinked_elastic/model.py',equationGroups:G1_V2_EQUATIONS,
  algorithm:['Generate 30 finite curved fibres around a 9 µm-radius cell and fix both ends of every fibre.','Detect geometric segment intersections between different fibres and create stress-free permanent links.','For each fibre, locate the closest segment point inside the 3 µm contact shell and 60° protrusion sector.','Normalize Gaussian weights only over eligible contacts and distribute the ramped total pull.','Compute stretching, reference-curvature bending, crosslink and cell-repulsion forces.','Advance free beads with the overdamped equation.','Report radial alignment, RMS displacement, energies, active power and drag power; compare linked and unlinked conditions.'],
  papers:[paper('lee','Supports finite bead–spring collagen fibres, fibre-scale lengths/thicknesses and inter-fibre linking.','The 30-fibre density, arbitrary force units or fixing every finite fibre end.'),paper('abhilash','Supports crosslink-connected, heterogeneous force transmission and radial fibre alignment around contractile cells.','That this small G1 network is percolated or quantitatively collagen-like.'),paper('wang','Supports tension-driven recruitment/alignment as a long-range transmission mechanism.','A purely visual change as sufficient validation of realignment.'),paper('notbohm','Motivates later tension/compression asymmetry and microbuckling tests.','The symmetric tensile/compressive springs actually used in G1 V2.')],
  results:[{name:'Direct test',value:'Same network and same force with the crosslink term switched off/on.'},{name:'Observable',value:'Near-cell radial order, RMS displacement and stored energies.'},{name:'Topology',value:'Permanent freely hinged links at segment intersections.'},{name:'Interpretation',value:'Qualitative transmission prototype; no physical-force calibration or percolation gate.'}],
  limitation:{heading:'HOWEVER · the network was over-anchored and visually ambiguous',text:'Fixing both ends of every finite fibre prevents the network from reorganizing as a connected material. The animation also hid most beads/crosslinks, so stretched lines looked like new protrusions rather than transmitted fibre motion.',items:['All fibre ends fixed, even when not near the outer boundary.','Only 30 fibres in a 100 µm field and no boundary-connectivity acceptance gate.','Symmetric tension/compression response misses fibre microbuckling.','Force unit is arbitrary.','First animation made beads, links and contact forces hard to identify.'],leadsTo:'G2 V2: 99 finite fibres, physical nN–µm calibration, boundary-only anchors, percolation gate, compression softening and explicit visual encodings.'},
  talk:{script:'“G1 V2 introduced the correct conceptual hierarchy: a hard contact region, Gaussian sharing inside it, and crosslinks as the only route for indirect fibres to feel the cell. The equation worked, but the boundary condition dominated: I fixed every fibre end. Therefore I treat it as a topology prototype, not the accepted network baseline.”',assumptions:['Crosslinks are permanent, stress-free and freely hinged.','One closest material contact per fibre.','Every finite fibre end represents the far field.'],nextChecks:['Anchor only beads that actually reach the outer boundary.','Require contact fibres to connect through the crosslink graph to the boundary.','Use physical EA and EI and tension/compression asymmetry.','Show every bead, link, contact and anchor explicitly.']}
};

NOTEBOOK_MODELS['g1-v3']={
  generation:'G1',version:'V3',status:'first fixed-versus-moving comparison',title:'Two-sided motor–clutch imbalance and cell translation',
  subtitle:'Opposing protrusions were added to ask whether stochastic left/right traction imbalance can move the cell while deforming the same collagen network.',
  question:'If the cell is released, does the reaction to unequal left/right clutch traction produce motion that is absent when the same cell center is fixed?',
  units:'Network geometry is in µm and time in seconds; clutch/network force remains an arbitrary G1 unit.',
  purpose:'比較同一個 collagen network、同一個 stochastic stream 下的 fixed cell 和 moving cell，讓 cell movement 只來自兩側 clutch traction 的不平衡。',
  brief:['「Two-sided clutch imbalance and cell movement；我要看 not moving and moving 兩種 condition。」','「stable imbalance pulling force from both sides changes the links and then the cell can move。」','「motor–clutch attachment sites 要跟著 fibres deform 重新更新。」'],
  chain:['left/right contact patches','stochastic clutch binding','traction imbalance','reaction force','cell translation'],
  simulation:{src:'versions/v3_two_sided_migration/demo/index.html',note:'The original comparison is preserved; its released-cell displacement is only about 0.18 µm and therefore difficult to see.'},
  source:REPO+'versions/v3_two_sided_migration/model.py',inherits:'g1-v2',
  equationGroups:[
    {title:'Two-sided stochastic motor–clutch module',equations:[
      {name:'Bound clutch spring force',tex:T`\[F_{s,i}=k_cx_{s,i},\qquad T_s=\sum_{i\in bound,s}F_{s,i}\]`,explain:'Each left/right effective clutch stores extension; side traction is the sum of bound clutch forces.'},
      {name:'Motor force–velocity relation',tex:T`\[v_{a,s}=v_0\max\!\left(0,1-\frac{T_s}{F_{stall}}\right)\]`,explain:'Actin retrograde flow slows as one side approaches stall force.',key:true},
      {name:'Relative clutch loading',tex:T`\[\dot x_{s,i}=\max(0,v_{a,s}-v_{sub,s})\]`,explain:'Clutch extension increases from actin motion relative to the contacted collagen material point.'},
      {name:'Bell rupture and binding probabilities',tex:T`\[k_{off,i}=k_{off}^0e^{F_i/F_b},\quad P_{off,i}=1-e^{-k_{off,i}\Delta t},\quad P_{on}=1-e^{-k_{on}\Delta t}\]`,explain:'Load accelerates unbinding; unbound clutches bind with a constant hazard.',key:true},
      {name:'Patch selection on binding',tex:T`\[\Pr(m\mid bind)=w_m\]`,explain:'A new clutch chooses one eligible material patch using the contact Gaussian weights.'},
      {name:'ECM reaction on the cell',tex:T`\[\mathbf R_c=-\sum_i\mathbf F_i^{act}\]`,explain:'Newton’s third-law reaction is the negative of the force placed on the ECM.'},
      {name:'G1 moving-cell law',tex:T`\[\zeta_c\dot{\mathbf r}_c=\mathbf R_c\]`,explain:'The released cell translates in full 2D; the fixed control sets its velocity to zero.',key:true},
      {name:'Cell position update',tex:T`\[\mathbf r_c^{n+1}=\mathbf r_c^n+\Delta t\,\mathbf R_c/\zeta_c\]`,explain:'No rotation, deformation, polarity memory or speed cap is included.'},
      {name:'Common indexed random stream',tex:T`\[U=H(seed,step,side,clutch,channel)\in[0,1)\]`,explain:'The fixed and moving conditions receive the same proposed random draws; mechanics can still change hazard outcomes.'}
    ]}
  ],
  algorithm:['Create identical crosslinked networks and two groups of 12 unbound clutches.','Find contact patches in left and right 60° sectors.','Use the indexed random stream to propose binding and rupture events.','Update actin velocity and clutch extensions from material-point motion relative to the cell.','Spread clutch forces through the same contact-patch interpolation onto collagen beads.','Compute the opposite reaction on the cell.','In the fixed condition set cell velocity to zero; in the released condition divide the full 2D reaction by cell drag.','Advance collagen beads and record cell center, clutch state, traction, alignment and displacement.'],
  papers:[paper('bell','Supports force-accelerated slip-bond off-rate.','A complete integrin catch/slip-bond description.'),paper('bangasser','Supports stochastic motor–clutch loading and motor force–velocity feedback.','The exact G1 parameter values or two-sided geometry.'),paper('prahl','Supports using opposing protrusion modules and recognizing that persistent polarity is needed for directional migration.','The unpolarized, six-second G1 trajectory as a migration calibration.'),paper('speed','Provides an order-of-magnitude migration-speed comparison for MDA-MB-231 in collagen.','The arbitrary G1 force/drag parameters or its 0.18 µm displacement.')],
  results:[{name:'Fixed condition',value:'Cell center remains exactly fixed by constraint.'},{name:'Released condition',value:'Only ~0.18 µm displacement in the published short trajectory, visually hard to distinguish.'},{name:'Controlled difference',value:'Same initial network and indexed random proposals; translation switch differs.'},{name:'Missing behavior',value:'No persistent polarity, speed calibration, axis restriction, rotation or crosslink breaking.'}],
  limitation:{heading:'HOWEVER · fixed and moving looked almost identical',text:'The six-second run and arbitrary cell drag produced too little cell displacement to explain visually. Full 2D translation also allowed incidental y-asymmetry to introduce a new mechanism before the intended left/right baseline was understood.',items:['Only ~0.18 µm motion against a 9 µm-radius cell.','No speed ensemble or experimental target.','No persistent side bias, so direction can average away.','Inherited over-anchored G1 network.','Full 2D translation confounds the two-protrusion question.'],leadsTo:'G2 V3: calibrated 30-minute trajectories, persistent polarity, one declared migration axis, common stream and a clearly visible 4.76 µm fixed/moving difference.'},
  talk:{script:'“G1 V3 established the reaction-force logic: clutch forces go into the ECM, and the opposite force divided by cell drag moves a released cell. The control was conceptually clean but visually weak and not calibrated. That failure told me to calibrate speed and constrain the first migration test to the left/right protrusion axis.”',assumptions:['Two protrusions are sufficient for the first migration decision.','The cell is rigid and translates without rotating or deforming.','G1 allows full 2D translation from any reaction component.'],nextChecks:['Introduce persistent polarity explicitly.','Calibrate path speed against a collagen migration range.','Run long enough for multi-micrometre motion.','Restrict the first released degree of freedom to the left/right axis.']}
};

NOTEBOOK_MODELS['g1-v4']={
  generation:'G1',version:'V4',status:'first plasticity hypothesis',title:'Stress-activated crosslink turnover and rest-state reset',
  subtitle:'A load–unload experiment tested whether breaking and reforming the same crosslink partners with a new rest state can retain deformation.',
  question:'Does force-accelerated turnover of selected crosslinks create residual alignment after the active cell force is removed?',
  units:'G1 geometry in µm/time in seconds with arbitrary force and kinetic scales.',
  purpose:'把 permanent remodeling 放最後：先在 elastic baseline 上加入 crosslink breaking/reforming，再用 unload 後的 residual 判斷是否真的永久。',
  brief:['「永久 remodeling 放最後；先不要在前面混入 crosslink breaking 的 complexity。」','「cell moves ultimately when links differ/break；先提出用什麼 theory 模擬。」','「比較 elastic control 和 plastic condition，force 離開後才判斷 remodeling。」'],
  chain:['fixed cell loads network','link force crosses threshold','dynamic link breaks','same partners rebind with new rest vector','unload residual'],
  simulation:{src:'versions/v4_plastic_remodeling/demo/index.html',note:'The original V4 animation compares permanent-elastic links with a same-partner break/re-form rule.'},
  source:REPO+'versions/v4_plastic_remodeling/model.py',inherits:'g1-v2',
  equationGroups:[
    {title:'G1 crosslink turnover rule',equations:[
      {name:'Bound-link force magnitude',tex:T`\[F_x=k_x\left\|(\mathbf p_b-\mathbf p_a)-\boldsymbol\ell_x^0\right\|\]`,explain:'Only the selected dynamic half of existing intersection links can turn over.'},
      {name:'Thresholded activation',tex:T`\[a_x=\frac{\max(0,F_x-F_{th})}{F_{scale}}\]`,explain:'Force below threshold does not add exponential activation.'},
      {name:'Force-accelerated break rate',tex:T`\[k_{off,x}=k_{off,x}^0e^{a_x},\qquad P_{break}=1-e^{-k_{off,x}\Delta t}\]`,explain:'Higher link force increases the breaking hazard.',key:true},
      {name:'Same-partner reformation',tex:T`\[P_{reform}=1-e^{-k_{on,x}\Delta t}\]`,explain:'An unbound link may reconnect the same two material neighborhoods.'},
      {name:'Rest-state reset',tex:T`\[\boldsymbol\ell_{x,new}^0=\mathbf p_b(t_{reform})-\mathbf p_a(t_{reform})\]`,explain:'Reformation erases the current link strain; this was the minimal plasticity mechanism.',key:true},
      {name:'Abrupt load/unload protocol',tex:T`\[F(t)=\begin{cases}F_{max}\min(1,t/t_{ramp}),&t<t_{load}\\0,&t\ge t_{load}\end{cases}\]`,explain:'The active force turns off abruptly after the loading interval.'},
      {name:'Residual alignment',tex:T`\[\Delta S_r^{res}=S_r(t_{final})-S_r(0)\]`,explain:'Only alignment retained after force removal is a remodeling candidate.'}
    ]}
  ],
  algorithm:['Generate the same G1 V2 network and mark 50% of original crosslinks as dynamic.','Run elastic and plastic conditions with the same initial network and loading.','At every step, compute the force on each bound dynamic link.','Use the thresholded exponential hazard to break links stochastically.','Allow an unbound link to re-form between the same partners and reset its rest vector to current separation.','Turn off the active force and continue overdamped relaxation.','Compare final alignment/displacement with the permanent-elastic control.'],
  papers:[paper('ban','Motivates testing crosslink-level mechanisms in load–unload plastic deformation of collagen networks.','That same-partner rest reset is the mechanism identified by Ban et al.'),paper('kim','Supports stress-induced irreversible rearrangement/sliding in dynamic collagen networks.','The specific threshold/off/on parameters used in G1 V4.'),paper('abhilash','Supports judging remodeling spatially and after load history.','That a nonzero finite-time residual automatically equals permanent plasticity.')],
  results:[{name:'Comparison',value:'Permanent-elastic versus stochastic same-partner rest-reset condition.'},{name:'Published outcome',value:'No resolved excess remodeling relative to the elastic control.'},{name:'Event counters',value:'Tracks bound state, total broken links and total reformed links.'},{name:'Interpretation',value:'Hypothesis test only; the rule does not create new fibre connectivity.'}],
  limitation:{heading:'HOWEVER · rest-resetting old links was the wrong topology test',text:'The rule reconnects the same two partners and can erase strain without representing new fibre–fibre contacts, sliding or merging. A finite-time residual also needs an elastic control before it can be called plastic.',items:['No new crosslink topology; same material neighborhoods rebind.','Half of original links are arbitrarily labeled dynamic.','Abrupt unloading and arbitrary kinetics.','Inherited G1 boundary and calibration flaws.','No gating on validated elastic transmission or migration.'],leadsTo:'G2 V4: gate plasticity after V2/V3, allow only genuinely new nearby/aligned contacts to form stress-free weak links, and compare excess residual against an elastic control.'},
  talk:{script:'“G1 V4 was deliberately a final-layer hypothesis. It asked whether force-triggered turnover plus a rest-state reset can retain deformation after unloading. The result was negative, and the mechanism was topologically weak because it reconnected the same partners. That motivated a stricter new-contact rule in G2.”',assumptions:['Exactly half of original links are dynamic.','Reformed links reconnect the same material locations.','A new rest vector is assigned instantly at reformation.'],nextChecks:['Require genuinely new inter-fibre contact and alignment.','Use gradual unload and an elastic finite-time control.','Do not claim plasticity unless residual exceeds the control.','Explore sliding/merging only as a separate future mechanism.']}
};

NOTEBOOK_MODELS['g2-v2']={
  generation:'G2',version:'V2',status:'validated elastic baseline',title:'Corrected crosslink force transmission',
  subtitle:'A physically scaled, boundary-connected 99-fibre network compares the same pull with and without crosslink forces while showing every bead, spring, link, contact and anchor.',
  question:'Under a fixed 5 nN pull, do permanent crosslinks increase intermediate-range response without numerical energy error, cell penetration or artificial end anchoring?',
  units:'Length µm, force nN, time s, modulus MPa. Cell radius 9 µm; fibre length 20–80 µm; fibre diameter 0.30 µm.',
  purpose:'修正 G1 的 over-anchoring 和不清楚的 visualization；讓其他 fibres 只能透過 crosslinks 受力，並用 near/middle/far metrics 驗證傳力。',
  brief:['「橘色點不要一直閃；我要知道它代表什麼。」','「所有 collagen 都要清楚看出是 bead-based：beads 連成 springs；crosslink 要明確畫出來。」','「被拉處不能像無端長出觸手；附近 collagen 應透過 crosslinks 有部分 movement/reorientation。」','「cell/collagen 尺度、finite fibre length、far boundary anchoring 都要合理。」'],
  chain:['boundary-connected finite fibres','contact shell + Gaussian weights','5 nN local pull','crosslink force transfer','shell displacement/alignment + energy check'],
  simulation:{src:'generations/g2_corrected/v2_crosslink_transmission/demo/index.html',note:'No-link and linked panels use the same 7,173-bead network. Force arrows are display-scaled shorter in this update; the underlying force values are unchanged.'},
  source:REPO+'generations/g2_corrected/common/model.py',equationGroups:G2_BASE_EQUATIONS,
  algorithm:['Generate 99 finite 20–80 µm fibres in a 180 µm field, including six near-cell contact fibres.','Discretize every fibre at 0.75 µm spacing and fix only beads inside the outer 2.5 µm boundary band.','Create freely hinged crosslinks at intersections and reject the network unless every contact fibre and ≥85% of all fibres connect to the boundary component.','Build two copies of the accepted network: one ignores crosslink force and one includes it.','Cache the closest eligible material contact on each fibre in the protrusion sector and normalize Gaussian weights.','Ramp the same total 5 nN pull, compute all elastic/repulsive forces and advance free beads overdamped.','Record beads, bond strain, link force, contact force, shell alignment/displacement, energy, power and minimum cell gap.','Validate convergence and compare linked versus no-link intermediate response.'],
  papers:[paper('lee','Supports fibre-scale geometry, wet-modulus scale and bead–spring/crosslink representation.','That one selected coarse-grained parameter set is a fit to a particular gel.'),paper('abhilash','Supports connected discrete-fibre alignment, heterogeneous deformation and long-range transmission.','That 21.6% is a universal biological effect size.'),paper('wang','Supports tension-driven fibre recruitment and alignment as a long-range mechanism.','That any visual line motion proves nonlinear recruitment without strain metrics.'),paper('notbohm','Motivates compression-softened fibres as a microbuckling candidate.','That collagen and fibrin have identical compression ratios.'),paper('han','Motivates near/intermediate/far spatial metrics around a contractile breast-cancer cell.','The exact 5 nN total pull or 2D geometry.')],
  results:[{name:'Accepted network',value:'99 fibres, 7,173 beads, 7,074 bonds, 383 crosslinks; all fibres boundary connected.'},{name:'Crosslink transmission',value:'Intermediate displacement is 21.6% larger with crosslink force than in the matched no-link control.'},{name:'Mechanical safety',value:'p99 bond strain ≈0.029%, maximum ≈0.104%, and minimum cell gap stays positive.'},{name:'Energy validation',value:'Relative overdamped energy-balance mismatch ≈9.2×10⁻⁵.'}],
  limitation:{heading:'Current limitation · elastic baseline, not complete collagen biology',text:'G2 V2 is the accepted mechanics baseline, but it still represents 2D finite fibres with linearized stretch/bend/link laws and a prescribed fixed-cell pull.',items:['No cell movement or stochastic clutches.','No crosslink breaking, sliding, merging or new-link formation.','No 3D entanglement, steric fibre–fibre contact or collagen viscoelastic calibration.','One parameter set is mechanism-first, not fitted to one experiment.'],leadsTo:'G2 V3 uses this exact validated ECM and adds two-sided stochastic clutches plus a controlled fixed/moving comparison.'},
  talk:{script:'“G2 V2 is my corrected elastic baseline. I use physical EA and EI, compression softening, boundary-only anchors and a percolation gate. The cell directly loads only contacts inside a 3 µm shell; Gaussian weights merely share a fixed 5 nN total. Other fibres can move only through crosslinks. I validate the response with shell displacement, strain and energy balance, not by making the animation dramatic.”',assumptions:['2D circular rigid cell and finite 2D fibres.','Permanent freely hinged links at initial intersections.','5 nN prescribed total pull and 1.5 µm Gaussian width.'],nextChecks:['Calibrate force/drag to one experimental gel and cell condition.','Test 3D network density and steric contacts.','Compare alignment and displacement profiles with imaging/force microscopy.']}
};

NOTEBOOK_MODELS['g2-v3']={
  generation:'G2',version:'V3',status:'validated migration baseline',title:'Corrected fixed versus moving cell',
  subtitle:'The validated G2 ECM is coupled to two 12-clutch protrusions; a persistent but stochastic side bias drives only the declared left/right translation degree of freedom.',
  question:'With the same ECM and indexed random stream, does releasing x-translation create visible, experimentally plausible cell motion while the fixed control stays exactly at zero?',
  units:'The inherited ECM uses µm–nN–s. Reported cell speed is converted to µm/min.',
  purpose:'把 G1 看不出來的 fixed/moving 差異修正到清楚但不誇張：使用 normal migration speed 範圍、30-min trajectory、common random stream 和一個 left/right axis。',
  brief:['「Fixed cell versus moving cell 看不出差異；在合理參數範圍內讓差別看得出來。」','「兩側 stable imbalance pulling force 導致 cell move；比較 moving 和 not moving。」','「attachment sites 要隨 cell 和 fibres 更新，不要固定在舊 beads。」'],
  chain:['persistent side-biased binding','Bell clutch turnover','left/right traction imbalance','ECM reaction','x-axis cell motion'],
  simulation:{src:'generations/g2_corrected/v3_two_sided_migration/demo/index.html',note:'Both panels share the corrected ECM and indexed random proposals. Thin clutch spokes replace large arrows in this view so the cell trajectory remains visible.'},
  source:REPO+'generations/g2_corrected/v3_two_sided_migration/model.py',inherits:'g2-v2',
  equationGroups:[
    {title:'Two-sided clutch loading and turnover',equations:[
      {name:'Bound clutch force and side traction',tex:T`\[F_{s,i}=b_{s,i}k_cx_{s,i},\qquad T_s=\sum_{i=1}^{N_c}F_{s,i}\]`,explain:'b is 1 when bound and 0 when unbound; N₍c₎=12 per side.'},
      {name:'Motor force–velocity relation',tex:T`\[v_{a,s}=v_0\max\!\left(0,1-\frac{T_s}{F_{stall}}\right)\]`,explain:'More clutch traction slows actin flow on that protrusion.',key:true},
      {name:'Weighted substrate speed at contacts',tex:T`\[v_{sub,s}=\sum_mw_m\left[(1-\alpha_m)\dot{\mathbf r}_{i_m}+\alpha_m\dot{\mathbf r}_{j_m}\right]\cdot\mathbf n_{in,m}\]`,explain:'The attachment environment follows the deforming material points rather than a fixed bead list.'},
      {name:'Clutch extension update',tex:T`\[x_{s,i}^{n+1}=x_{s,i}^{n}+\Delta t\max(0,v_{a,s}-v_{sub,s})\quad\text{if bound}\]`,explain:'Relative actin/substrate motion loads a bound clutch.'},
      {name:'Bell slip-bond off-rate',tex:T`\[k_{off,s,i}=k_{off}^0\exp(F_{s,i}/F_b),\qquad P_{off}=1-e^{-k_{off,s,i}\Delta t}\]`,explain:'High force shortens a clutch lifetime.',key:true},
      {name:'Persistent polarity in binding hazard',tex:T`\[k_{on,L}=2(1-\psi)k_{on},\qquad k_{on,R}=2\psi k_{on},\qquad \psi=0.65\]`,explain:'The mean on-rate is conserved while a persistent side preference prevents purely zero-mean wandering.',key:true},
      {name:'Binding probability',tex:T`\[P_{on,s}=1-e^{-k_{on,s}\Delta t}\]`,explain:'The common indexed random array is compared with side-specific probabilities.'},
      {name:'Clutch force distributed to ECM',tex:T`\[\mathbf F^{act}_{m,s}=T_sw_m\mathbf n_{in,m}\]`,explain:'Side traction is distributed through the current contact patches with the same G2 contact weights.'}
    ]},
    {title:'Cell reaction, constrained translation and speed calibration',equations:[
      {name:'Reaction force',tex:T`\[\mathbf R_c=-\sum_i\mathbf F_i^{act}\]`,explain:'The ECM force has an equal and opposite reaction on the cell.'},
      {name:'One-axis cell law',tex:T`\[\dot x_c=\operatorname{clip}\!\left(\frac{R_{c,x}}{\zeta_c},-v_{max},v_{max}\right),\qquad \dot y_c=0\]`,explain:'V3 intentionally releases only the left/right protrusion axis; the fixed control sets both components to zero.',key:true},
      {name:'Cell-center update',tex:T`\[x_c^{n+1}=x_c^n+\Delta t\,\dot x_c\]`,explain:'The contact patches are periodically recomputed around the new center.'},
      {name:'Common random stream',tex:T`\[U^{on/off}_{n,s,i}=H(seed,n,s,i,channel)\]`,explain:'Fixed and moving conditions receive identical proposed random numbers; only their mechanics can make the trajectories diverge.'},
      {name:'Path/instantaneous speed calibration metric',tex:T`\[\bar v_{path}=\frac{60}{T}\sum_n|\dot x_c^n|\Delta t\quad[\mu m/min]\]`,explain:'The ensemble target is median 0.2–0.4 µm/min, comparable to MDA-MB-231 collagen measurements.',key:true},
      {name:'Net displacement and instantaneous speed',tex:T`\[\Delta x_c=x_c(T)-x_c(0),\qquad v_{inst}=60|\dot x_c|\]`,explain:'The animation reports both accumulated position and current speed.'}
    ]}
  ],
  algorithm:['Load one accepted G2 V2 network and clone it for fixed and moving conditions.','Generate the same indexed on/off random arrays for both conditions.','Initialize 12 clutches on each side and find left/right contact patches.','Compute motor speed, substrate speed, extension, Bell rupture and side-biased binding.','Distribute left and right side traction onto current material contact points.','Apply the opposite reaction to the cell; constrain the fixed center or update only x for the released cell with a speed cap.','Advance the same crosslinked ECM mechanics and periodically update contact patches as fibres/cell move.','Record 6 s frames for 30 min and separately run a 20-seed speed ensemble.'],
  papers:[paper('bell','Supports the exponential force-dependent slip-bond off-rate.','Integrin catch-bond complexity or side polarity.'),paper('bangasser','Supports stochastic clutch loading and motor force–velocity feedback.','The exact 12-clutch and drag parameters.'),paper('prahl','Supports opposing protrusion modules and persistent polarity for directional migration.','That the present rigid-cell one-axis model is a complete migration model.'),paper('speed','Provides an MDA-MB-231 speed scale in collagen that brackets the 0.2–0.4 µm/min calibration target.','A fit to one collagen concentration, cell line state or trajectory.'),paper('steinwachs','Supports using breast-cancer-cell forces in collagen as a quantitative future calibration context.','The current cell drag value without direct fitting.')],
  results:[{name:'20-seed speed ensemble',value:'Median path speed 0.280 µm/min; seed range about 0.251–0.312 µm/min.'},{name:'Fixed cell',value:'Exactly 0 µm cell displacement by constraint.'},{name:'Moving cell',value:'4.76 µm net x displacement and 8.18 µm path over 30 min.'},{name:'Mechanical safety',value:'Global maximum bond strain remains below the 15% validation gate and the cell–fibre gap stays positive.'}],
  limitation:{heading:'Current limitation · directional migration is constrained, not emergent in 2D',text:'The model demonstrates how a declared persistent polarity plus clutch traction can move a rigid cell along one axis. It does not yet explain how polarity itself emerges, nor does it include rotation or deformation.',items:['ψ=0.65 is prescribed.','Only x translation is released; y, rotation and shape are fixed.','Crosslinks remain permanent.','Cell drag and clutch parameters are mechanism-calibrated, not jointly fitted.'],leadsTo:'G2 V4 keeps the cell fixed again and isolates a separate question: can load create new inter-fibre contacts that remain after unloading?'},
  talk:{script:'“G2 V3 changes only the cell constraint between the two panels. Both see the same ECM and random proposals. The clutches load according to actin minus substrate speed, rupture by Bell kinetics and exert the opposite reaction on the cell. I prescribe a modest right-side binding bias and release only x, which gives a median 0.28 µm/min and a visible 4.76 µm displacement without inventing 2D polarity.”',assumptions:['Rigid circular cell with one translational degree of freedom.','Persistent polarity ψ=0.65 is prescribed rather than generated.','Slip bonds, 12 effective clutches per side and a speed cap.'],nextChecks:['Fit clutch and cell-drag parameters jointly to force and speed data.','Add rotation and y translation only after the one-axis result is robust.','Model polarity dynamics rather than prescribing ψ.','Later couple migration to validated crosslink turnover/breaking.']}
};

NOTEBOOK_MODELS['g2-v4']={
  generation:'G2',version:'V4',status:'gated plasticity candidate · negative baseline',title:'New-contact weak links under load–unload',
  subtitle:'A stricter plasticity candidate allows only genuinely approaching, aligned fibres to form new stress-free weak links and compares the fully unloaded residual with an elastic control.',
  question:'Does a validated 5 nN elastic loading history create enough genuinely new fibre contact to leave excess residual deformation after complete unloading?',
  units:'Inherited G2 µm–nN–s mechanics; 300 s load, 30 s unload ramp and 900 s total observation.',
  purpose:'不要把 loading deformation 當 permanent remodeling；只有 force 完全移除後、而且超過 elastic control 的 residual 才算候選 plasticity。',
  brief:['「永久 remodeling 放最後；crosslink breaking/new topology 是另一層 complexity。」','「先讓 elastic baseline 完整，再測 force 離開後怎麼辦。」','「不能為了動畫明顯就任意調大參數；negative result 也要誠實顯示。」'],
  chain:['validated elastic load','different fibres approach + align','new stress-free weak link','gradual unload','excess residual versus elastic control'],
  simulation:{src:'generations/g2_corrected/v4_contact_plasticity/demo/index.html',note:'The left panel is elastic; the right can form new purple weak links. Contact arrows have been shortened; their numerical forces are unchanged.'},
  source:REPO+'generations/g2_corrected/v4_contact_plasticity/model.py',inherits:'g2-v2',
  equationGroups:[
    {title:'Load–unload protocol',equations:[
      {name:'Piecewise force history',tex:T`\[F(t)=\begin{cases}F_{max}t/t_r,&0\le t<t_r\\F_{max},&t_r\le t<t_L\\F_{max}\left[1-\frac{t-t_L}{t_u}\right],&t_L\le t<t_L+t_u\\0,&t\ge t_L+t_u\end{cases}\]`,explain:'The force ramps up, holds, ramps down over 30 s and remains fully zero for the final interval.',key:true},
      {name:'Elastic control',tex:T`\[\mathcal M_{elastic}:\quad N_{new}(t)=0\]`,explain:'The control uses the same network and force protocol but cannot form new links.'}
    ]},
    {title:'New-contact formation criteria',equations:[
      {name:'Different-fibre proximity',tex:T`\[f_a\ne f_b,\qquad d_{ab}(t)=\|\mathbf r_b-\mathbf r_a\|\le d_{cap}=0.45\ \mu m\]`,explain:'Only beads belonging to different fibres can form a candidate link.'},
      {name:'Approach under load',tex:T`\[d_{ab}(0)-d_{ab}(t)\ge\Delta d_{min}=2.5\times10^{-4}\ \mu m\]`,explain:'The pair must become closer than it was initially; pre-existing near neighbors do not count.',key:true},
      {name:'Alignment requirement',tex:T`\[|\hat{\mathbf t}_a\cdot\hat{\mathbf t}_b|\ge\cos30^\circ\]`,explain:'The fibres must be approximately parallel, allowing either tangent orientation.'},
      {name:'Near-cell formation region',tex:T`\[0\le\|\mathbf r_{a,b}-\mathbf r_c\|-R\le30\ \mu m\]`,explain:'New links are considered only in the mechanically affected neighborhood.'},
      {name:'No duplicate local topology',tex:T`\[\|\mathbf m_{candidate}-\mathbf m_{existing}\|\ge1\ \mu m\quad\text{for the same fibre pair}\]`,explain:'This prevents stacking several links at the same local contact.'},
      {name:'Stress-free new rest vector',tex:T`\[\boldsymbol\ell_{new}^0=\mathbf p_b(t_f)-\mathbf p_a(t_f),\qquad \mathbf F_{new}(t_f)=0\]`,explain:'A newly formed weak link stores no force at the deformed geometry where it forms.',key:true},
      {name:'Weak-link mechanics',tex:T`\[\mathbf F_{new}=k_{weak}\left[(\mathbf p_b-\mathbf p_a)-\boldsymbol\ell_{new}^0\right],\qquad k_{weak}=15\ \mathrm{nN/\mu m}\]`,explain:'New links use the same material-point interpolation as permanent links but lower stiffness.'},
      {name:'Formation-rate caps',tex:T`\[N_{new/check}\le3,\qquad N_{new,total}\le36,\qquad \Delta t_{check}=5\ s\]`,explain:'Deterministic caps keep one geometric check from creating a dense bundle instantly.'}
    ]},
    {title:'Plasticity decision metrics',equations:[
      {name:'Finite-time unloaded residual',tex:T`\[u_{res}(t)=\sqrt{N_{free}^{-1}\sum_i\|\mathbf r_i(t)-\mathbf r_i^0\|^2},\qquad t\ge t_L+t_u\]`,explain:'This is measured only after the active force reaches zero.'},
      {name:'Excess residual over elastic control',tex:T`\[\Delta u_{plastic}=u_{res}^{new-link}-u_{res}^{elastic}\]`,explain:'Plasticity requires an excess above the elastic finite-time residual, not merely a nonzero displacement.',key:true},
      {name:'Excess alignment residual',tex:T`\[\Delta S_{plastic}=(S_r^{new-link}(t_f)-S_r(0))-(S_r^{elastic}(t_f)-S_r(0))\]`,explain:'The same control subtraction applies to alignment.'},
      {name:'Stage gate',tex:T`\[G_{V4}=G_{V2,network}\land G_{V2,mechanics}\land G_{V3,mobility}\]`,explain:'V4 refuses to run unless the corrected network, mechanics and mobility validations have passed.'}
    ]}
  ],
  algorithm:['Check that G2 V2 network/mechanics and G2 V3 mobility gates are true.','Create identical elastic and new-link networks from one accepted G2 specification.','Apply the 5 nN ramp–hold–ramp-down force protocol and update current contact patches.','Every 5 s during loading, search nearby nonfixed beads on different fibres.','Require new approach, ≤0.45 µm distance, ≤30° relative alignment, near-cell location and no duplicate local link.','Create at most three stress-free 15 nN/µm weak links per check; do not break or reset original links.','Continue the overdamped simulation through a long fully unloaded interval.','Compare residual RMS and alignment against the elastic control before making any plasticity claim.'],
  papers:[paper('ban','Motivates load-dependent formation of new weak inter-fibre links as a candidate collagen plasticity mechanism.','That this minimal deterministic capture rule reproduces the full model or experiment.'),paper('kim','Supports irreversible sliding/merging as an alternative dynamic-network mechanism.','That sliding/merging is implemented here—it is not.'),paper('abhilash','Supports evaluating retained spatial remodeling rather than only local loading deformation.','That a finite-time residual without a control is permanent.'),paper('lee','Supports material-point fibre/crosslink representation and geometric scales.','The selected weak-link capture and stiffness parameters.')],
  results:[{name:'New links formed',value:'One genuinely new weak link in the accepted baseline trajectory.'},{name:'Elastic unloaded RMS',value:'0.0051287 µm at the final observation time.'},{name:'New-link unloaded RMS',value:'0.0051287 µm at the final observation time.'},{name:'Conclusion',value:'No resolved excess residual; this baseline does not support a permanent-remodeling claim.'}],
  limitation:{heading:'Current limitation · honest negative result',text:'At the accepted 5 nN load, the elastic network rarely brings new fibre pairs together. One new weak link is not enough to create measurable excess residual. The visualization therefore remains subtle.',items:['No original-link breaking or sliding.','Capture/stiffness thresholds are mechanism-first assumptions.','Only one loading history and 2D geometry.','Finite observation time, although controlled by an elastic comparison.'],leadsTo:'A future G3 should vary loading duration/history or implement a separately justified sliding/merging/breaking mechanism—without silently tuning the animation.'},
  talk:{script:'“G2 V4 is a gated negative experiment. A pair must be on different fibres, become newly closer, be within 0.45 µm and aligned within 30°. The new link is stress-free when created. After gradual unload I subtract the elastic residual. The baseline forms one link but shows no excess residual, so I do not claim permanent remodeling.”',assumptions:['New-link formation is deterministic once geometric criteria are met.','Original permanent links cannot break or reset in this stage.','A stress-free weak spring is an adequate first candidate for a new contact.'],nextChecks:['Sweep only justified loading history and plastic mechanism parameters.','Compare new-link formation frequency with imaging if available.','Implement sliding/merging or breakage as separate hypotheses with their own controls.','Extend the validated load–unload test to 3D.']}
};
