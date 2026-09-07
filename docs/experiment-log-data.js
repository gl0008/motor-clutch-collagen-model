window.REPOSITORY_BRANCHES = [
  {name:'main',role:'Complete synthesis',state:'canonical',description:'Contains the reviewed code, documentation and public research notebook. Merging here does not delete any generation branch.'},
  {name:'generation/g1',role:'Generation 1 checkpoint',state:'preserved',description:'Permanent G1 branch. The G1 release tags identify its individual historical versions.'},
  {name:'generation/g2',role:'Generation 2 checkpoint',state:'preserved',description:'Permanent corrected G2 branch, kept separately after integration into main.'},
  {name:'generation/g3',role:'Generation 3 checkpoint',state:'preserved',description:'Permanent G3 spheroid-guidance branch at the accepted G3 evidence state.'},
  {name:'generation/g4',role:'Generation 4 checkpoint',state:'preserved',description:'Permanent G4 branch at the long-time and multiscale v2 evidence state.'},
  {name:'generation/g5',role:'Generation 5 checkpoint',state:'preserved',description:'Permanent G5 model branch at the current molecular-clutch and ablation evidence state.'},
  {name:'agent/g5-organoid-plan',role:'G5 working history',state:'active history',description:'The model-development branch remains visible after merge so its experiment sequence can be audited.'},
  {name:'codex/model-evolution-notebook',role:'Notebook working history',state:'active history',description:'The website-development branch remains visible after merge so documentation changes can be audited.'}
];

window.VERSION_EVIDENCE = {
  'g1-v0':{branch:'generation/g1',commit:'e637d47'},'g1-v1':{branch:'generation/g1',commit:'e637d47'},
  'g1-v2':{branch:'generation/g1',commit:'e637d47'},'g1-v3':{branch:'generation/g1',commit:'e637d47'},
  'g1-v4':{branch:'generation/g1',commit:'e637d47'},
  'g2-v2':{branch:'generation/g2',commit:'96f6cc1'},'g2-v3':{branch:'generation/g2',commit:'96f6cc1'},
  'g2-v4':{branch:'generation/g2',commit:'96f6cc1'},
  'g3':{branch:'generation/g3',commit:'c795bfd'},
  'g4-v1-a':{branch:'agent/g4-interactive-calibration',commit:'5922afa'},
  'g4-v1-b':{branch:'agent/g4-interactive-calibration',commit:'5922afa'},
  'g4-v1-c':{branch:'agent/g4-interactive-calibration',commit:'5922afa'},
  'g4-v1-d':{branch:'agent/g4-interactive-calibration',commit:'5922afa'},
  'g4-v2-a':{branch:'generation/g4',commit:'de0c0ba'},'g4-v2-b':{branch:'generation/g4',commit:'de0c0ba'},
  'g4-v2-c':{branch:'generation/g4',commit:'de0c0ba'},'g4-v2-d':{branch:'generation/g4',commit:'de0c0ba'},
  'g5-0a':{branch:'generation/g5',commit:'2227642'},'g5-0b':{branch:'generation/g5',commit:'8287490'},
  'g5-0c':{branch:'generation/g5',commit:'60f7a75'},'g5-a':{branch:'generation/g5',commit:'3d2706a'},
  'g5-b':{branch:'generation/g5',commit:'3d2706a'},'g5-c':{branch:'generation/g5',commit:'9f5f68b'},
  'g5-d':{branch:'generation/g5',commit:'d6bce90'},'g5-d2':{branch:'generation/g5',commit:'1651038'},
  'g5-d3':{branch:'generation/g5',commit:'a6b3ae2'},'g5-e':{branch:'generation/g5',commit:'7d2bda9'},
  'g5-ablation':{branch:'generation/g5',commit:'60f7a75'}
};

window.EXPERIMENT_LOG = [
  {
    date:'2026-09-07',generation:'G3–G5',version:'Theory record',branch:'codex/model-evolution-notebook',commit:'c4085e0',
    parent:'Validated implementations and result pages',
    question:'Can the implemented mechanics be reconstructed without reading the source code?',
    change:'Typeset separate English theory guides for G3, G4 and G5.',
    result:'Each generation now has an equation-level companion document linked from its evidence page.',
    limitation:'The guides document the model; they do not add validation data.',
    next:'Keep every future theory guide synchronized with its accepted generation checkpoint.',
    evidence:[{label:'G3 theory',href:'theory/g3_model_theory.pdf'},{label:'G4 theory',href:'theory/g4_model_theory.pdf'},{label:'G5 theory',href:'theory/g5_model_theory.pdf'}]
  },
  {
    date:'2026-09-04',generation:'G5',version:'Ablation evidence set',branch:'agent/g5-organoid-plan',commit:'bba8f40',
    parent:'G5 one-factor ablation',
    question:'Can each mechanism be inspected without shrinking eight panels into one figure?',
    change:'Split the combined ablation result into eight mechanism-specific animations and retained the combined comparison.',
    result:'Baseline, no-clutch, no-adhesion, no-crosslink, stiffening, plasticity, leader and high-adhesion evidence are independently inspectable.',
    limitation:'The runs still use one two-dimensional seed and parameter set.',
    next:'Repeat the most informative mechanisms across seeds and experimental calibration targets.',
    evidence:[{label:'Combined ablation',href:'g5.html#clutch-baseline'}]
  },
  {
    date:'2026-09-04',generation:'G5',version:'G5-0C + ablation',branch:'agent/g5-organoid-plan',commit:'60f7a75',
    parent:'G5-0B cell-number control',
    question:'Does adhesion generate motion, or does it keep clutch-driven cells together?',
    change:'Released matched cells with adhesion off and compared them with the adhesion-on molecular-clutch baseline; also completed the eight-condition ablation.',
    result:'Adhesion-off cells move 68 µm and disperse; adhesion-on cells move 25 µm as a cohesive front. Removing the clutch or crosslinks sharply reduces motion.',
    limitation:'The result identifies internal causal roles but is not yet calibrated to an organoid movie.',
    next:'Calibrate clutch force, adhesion and collagen transmission against measured motion and cohesion.',
    evidence:[{label:'G5-0C and ablation',href:'g5.html#clutch-baseline'}]
  },
  {
    date:'2026-09-04',generation:'G5',version:'G5-0B',branch:'agent/g5-organoid-plan',commit:'8287490',
    parent:'G5-0A scale transfer',
    question:'Does increasing cell number amplify the traction of an individual active clutch sector?',
    change:'Varied only cell number on a corona collagen network that restored angular fibre access.',
    result:'Per-active-sector traction remains 4.40–4.51 nN from 1 to 19 cells; there is no per-clutch collective amplification.',
    limitation:'Fixed cells do not reveal whether adhesion changes motility or only cohesion.',
    next:'Release the cells with adhesion switched off while keeping the clutch unchanged.',
    evidence:[{label:'Lineage controls',href:'g5.html#lineage-controls'}]
  },
  {
    date:'2026-09-04',generation:'G5',version:'G5-0A',branch:'agent/g5-organoid-plan',commit:'2227642',
    parent:'G4 single-cell clutch',
    question:'Does G4 clutch physics transfer to the larger G5 domain when non-multicellular parameters are frozen?',
    change:'Ran one cell in the G5 box with the G4 parameter set and matched fibre density.',
    result:'Per-active-sector traction transfers within 15%; lower total traction is caused by sparse contact sectors.',
    limitation:'One cell cannot test whether cell number or adhesion produces collective behaviour.',
    next:'Keep the frozen physics and vary cell number only.',
    evidence:[{label:'G5-0A control',href:'g5.html#lineage-controls'}]
  },
  {
    date:'2026-09-04',generation:'G5',version:'G5D v3',branch:'agent/g5-organoid-plan',commit:'a6b3ae2',
    parent:'G5D v2 low-adhesion leaders',
    question:'Does extra leader traction pull followers into a connected strand?',
    change:'Added leader traction as a separate factor instead of bundling it with reduced adhesion.',
    result:'The leader remains fast, but followers do not advance as a connected strand in the tested regime.',
    limitation:'The experiment used the earlier prescribed-pull baseline.',
    next:'Preserve the negative result and retest mechanisms around the molecular-clutch baseline.',
    evidence:[{label:'Earlier leader test',href:'g5.html#leader-update'}]
  },
  {
    date:'2026-09-03',generation:'G5',version:'G5D v2',branch:'agent/g5-organoid-plan',commit:'1651038',
    parent:'G5D homogeneous invasion',
    question:'Can reduced adhesion expose cohesive, collective and single-cell invasion modes?',
    change:'Assigned a sparse outer leader population reduced adhesion and mapped adhesion against pull.',
    result:'The model separates connected motion from single-cell escape, but the apparent leader often outruns its followers.',
    limitation:'Adhesion and prescribed traction were still mixed across the original mode map.',
    next:'Separate leader traction from leader adhesion.',
    evidence:[{label:'Mode history',href:'g5.html#leader-update'}]
  },
  {
    date:'2026-09-03',generation:'G5',version:'Molecular-clutch baseline',branch:'agent/g5-organoid-plan',commit:'d6bce90',
    parent:'G5 A–E constant-pull implementation',
    question:'Can G5 inherit G4 Bell-law clutch loading and failure instead of prescribing constant traction?',
    change:'Added cell–fibre molecular-clutch dynamics, validation and explicit clutch reporting.',
    result:'G5 traction becomes emergent through load-and-fail events and can be used as the biological baseline.',
    limitation:'Changing the force mechanism requires lineage controls before comparing G4 and G5.',
    next:'Freeze the G4 parameters and audit scale transfer before attributing changes to multicellularity.',
    evidence:[{label:'Current G5 baseline',href:'g5.html#clutch-baseline'}]
  },
  {
    date:'2026-09-03',generation:'G5',version:'G5 A–E',branch:'agent/g5-organoid-plan',commit:'3d2706a',
    parent:'G4 reaction-driven single cell',
    question:'What changes when many adhesive cells pull one shared collagen network?',
    change:'Integrated the organoid scaffold, collective contraction, stiffness test, released-cell invasion and crosslink plasticity into one staged model.',
    result:'The complete multicellular implementation runs efficiently and exposes cohesion, invasion and remodelling readouts.',
    limitation:'Several parameters changed together relative to G4, so the raw generational difference was not attributable.',
    next:'Replace the constant-pull default with the molecular clutch and add explicit lineage controls.',
    evidence:[{label:'G5 staged model',href:'g5.html'}]
  },
  {
    date:'2026-09-01',generation:'G4',version:'G4 v2 A–D',branch:'generation/g4',commit:'de0c0ba',
    parent:'G4 v1 short-time laboratory',
    question:'Can longer observation and multiscale views reveal physical movement without changing the geometry?',
    change:'Extended A, B and D to two hours; added 1×, 10× and 50× displacement views and a clutch-event microscope.',
    result:'G4 now separates long-horizon ECM response, graph-path transmission, clutch failure and reaction-driven motion.',
    limitation:'A single rigid 2D cell cannot represent organoid cohesion or collective invasion.',
    next:'Create a separate multicellular generation rather than enlarging G4 silently.',
    evidence:[{label:'G4 v2 laboratory',href:'g4-v2.html'},{label:'G4 theory',href:'theory/g4_model_theory.pdf'}]
  },
  {
    date:'2026-08-24',generation:'G4',version:'G4 v1 A–D',branch:'agent/g4-interactive-calibration',commit:'5922afa',
    parent:'G3 radial-remodelling demonstration',
    question:'Which ECM and clutch mechanisms independently control transmission, slippage and movement?',
    change:'Separated mechanics calibration, crosslink-path transmission, Bell slippage and released-cell motion into four controlled stages.',
    result:'The model obtained a clean experimental sequence, but short observation made small motions difficult to read.',
    limitation:'Two-to-five-minute views hid the long-time response and individual failure events.',
    next:'Preserve v1 and build a longer multiscale G4 v2.',
    evidence:[{label:'Preserved G4 v1',href:'g4-lab.html'}]
  },
  {
    date:'2026-08-18',generation:'G3',version:'Spheroid guidance',branch:'generation/g3',commit:'c795bfd',
    parent:'G2 prescribed-polarity migration',
    question:'Can unbiased all-around protrusions reorganise collagen radially without choosing a migration direction?',
    change:'Removed prescribed polarity, introduced a spheroid gap and symmetric grip-and-reel protrusions, and aligned the laboratory with conservative contact.',
    result:'The near-field network moves toward a radial aster while the spheroid remains nearly fixed.',
    limitation:'The demonstration uses softened collagen and does not isolate transmission or clutch failure.',
    next:'Calibrate the mechanisms one at a time in G4.',
    evidence:[{label:'G3 model notebook',href:'model-notebook.html?model=g3'},{label:'G3 theory',href:'theory/g3_model_theory.pdf'}]
  },
  {
    date:'2026-08-14',generation:'G2',version:'Corrected G2 baseline',branch:'generation/g2',commit:'96f6cc1',
    parent:'G1 question-forming prototypes',
    question:'Can the original questions be rebuilt with physical units, outer-boundary anchoring and measurable controls?',
    change:'Built corrected crosslink transmission, calibrated constrained migration and a gated new-contact plasticity control.',
    result:'The network passes transmission and mobility checks; the plasticity control remains an honest negative result.',
    limitation:'Migration still used prescribed polarity, and the network response remained too stiff and local for the organoid movies.',
    next:'Remove preset direction and let spatial protrusions generate radial organisation.',
    evidence:[{label:'G2 version library',href:'index.html#versions'}]
  },
  {
    date:'2026-08-10',generation:'G1',version:'G1 V0–V4',branch:'generation/g1',commit:'e637d47',
    parent:'Initial model question',
    question:'Which minimal elements are needed to connect cell traction, collagen transmission, migration and plasticity?',
    change:'Created the SLS, few-fibre, crosslinked, moving-cell and first plasticity prototypes as separate preserved versions.',
    result:'The prototypes established the project questions and exposed ambiguous units, boundaries and visual interpretation.',
    limitation:'Prototype assumptions were added before a corrected collagen baseline was validated.',
    next:'Preserve G1 and rebuild its questions as Generation 2.',
    evidence:[{label:'G1 version library',href:'index.html#versions'}]
  }
];
