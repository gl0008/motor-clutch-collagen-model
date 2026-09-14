# G4 v4 assumption ledger

| ID | Assumption | Why required | Expected impact | Check |
|---|---|---|---|---|
| V4-A001 | The primary system is one MDA-MB-231 cell. | It gives a source-matched two-hour displacement and six-hour migration clock. | Results do not calibrate CT26 spheroids or other cell lines. | Keep spheroid work in a separate future configuration. |
| V4-A002 | Cell radius is provisionally 10 µm and constant. | No segmented cell geometry was supplied. | Contact number and force leverage depend on radius. | Repeat at 8, 10, and 12 µm and report distances normalized by radius. |
| V4-A003 | The model remains 2D. | It isolates contact, rotation, stretching, bending and crosslink paths. | Geometric intersections and connectivity are more frequent than in 3D. | Treat mechanism direction, not fitted crosslink probability, as transferable. |
| V4-A004 | Initial source-fiber angles are uniform on \([0,\pi)\). | Alignment must be an output. | Finite shells still have nonzero initial order. | Report \(S_r(0)\) and use \(\Delta S_r\). |
| V4-A005 | Fibers start straight and stress free. | It provides a controlled bending-to-stretching test. | Missing crimp may recruit axial tension too early. | Add measured waviness only after an experiment is selected. |
| V4-A006 | Crosslinks are permanent elastic material-point links. | The revision isolates elastic force transmission. | It cannot create permanent remodeling after force removal. | Crosslink breaking/plasticity remain outside G4 v4. |
| V4-A007 | A 2D intersection is retained with probability \(p_x\). | It creates nested connectivity without changing geometry. | \(p_x\) is not a biochemical crosslinker concentration. | Report candidate/retained links and graph connectivity. |
| V4-A008 | Outer-boundary beads are fixed. | It prevents rigid translation of the finite ECM. | A nearby boundary can stiffen the response. | Compare increasing domains while preserving network density. |
| V4-A009 | Direct force uses a 0–3 µm band and \(\sigma_c=1.5\) µm Gaussian weights. | These are the advisor-approved inherited contact rules. | Contact number depends on local random geometry. | Keep values visibly provisional until a matching experiment is found. |
| V4-A010 | Normal surface traction is primary. | It follows the professor’s proposed minimal model. | It is less anisotropic than an elongated cell. | Compare a two-sector dipole without changing total force. |
| V4-A011 | The dipole sector half-width is 30°. | A finite sector is required on a discrete random network. | The exact contact count varies by seed. | Renormalize force and report retained contacts. |
| V4-A012 | The primary site uses equal load sharing; independent Bell clutches remain a control. | Independent redundancy produced many molecular ruptures but no complete site failures in the pipeline check. | Equal sharing raises remaining-clutch load as the bound count falls and can create a cascade. Integrin catch bonding, reinforcement and maturation remain omitted. | Show both cases with identical ECM and random stream; label constants “effective.” |
| V4-A013 | A contact material point changes only after every clutch at that site is unbound. | It prevents attachment teleportation while load is transmitted. | A fully detached site may relocate discontinuously. | Log every site failure and relocation with the same timestamp. |
| V4-A014 | After complete site failure, the outgoing direction is the encountered local fiber tangent oriented away from the cell; the eligible surface patch nearest that direction is selected. | Riching et al. support collagen-guided protrusion orientation, but do not provide a relocation probability law. | It can bias successive contacts along locally encountered collagen without a global polarity vector. | Treat the selection rule as an explicit assumption and compare it with a no-guidance control later. |
| V4-A015 | Cell translation is overdamped; the circular cell has no rotational state. | Rotation of an isotropic circle is unobservable. | Shape-driven torque and deformation are absent. | Add rotation only with an anisotropic or deformable cell. |
| V4-A016 | New contacts exist only where collagen is encountered. | It avoids prescribed global polarity. | Geometry and stochastic lifetime can create force imbalance. | Compare fixed/released runs with identical seed and contact stream. |
| V4-A017 | Brownian force, SLS, plasticity, crosslink rupture, growth and degradation are absent. | Each is a separate mechanism that would obscure this elastic/contact test. | Long-term residual remodeling is not represented. | List them in every limitation/“However” section. |

## Parameter status

The two-hour and six-hour observation clocks are source matched. The 1 and
4 mg/ml collagen conditions in Riching et al. are documented comparison
conditions, but this bead network is not yet quantitatively calibrated to either
gel. The inherited stiffness, drag, force, contact and clutch values remain
mechanism-first defaults.
