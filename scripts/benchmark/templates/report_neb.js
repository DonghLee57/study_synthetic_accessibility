    (function initNebTrajectoryExplorers() {
      const CPK_STYLE = {
        Zr: { fill: '#0891b2', hi: '#67e8f9', rim: '#0e7490', cov: 1.75, rad: 17.5, txt: '#ffffff' },
        Ti: { fill: '#4f46e5', hi: '#a5b4fc', rim: '#3730a3', cov: 1.60, rad: 16.5, txt: '#ffffff' },
        Ga: { fill: '#c026d3', hi: '#f5d0fe', rim: '#86198f', cov: 1.22, rad: 15.5, txt: '#ffffff' },
        Si: { fill: '#d97706', hi: '#fde68a', rim: '#92400e', cov: 1.11, rad: 14.5, txt: '#ffffff' },
        Cl: { fill: '#16a34a', hi: '#86efac', rim: '#166534', cov: 1.02, rad: 12.8, txt: '#ffffff' },
        O:  { fill: '#dc2626', hi: '#fca5a5', rim: '#991b1b', cov: 0.66, rad: 11.2, txt: '#ffffff' },
        N:  { fill: '#2563eb', hi: '#93c5fd', rim: '#1e40af', cov: 0.71, rad: 11.6, txt: '#ffffff' },
        C:  { fill: '#334155', hi: '#94a3b8', rim: '#0f172a', cov: 0.76, rad: 10.8, txt: '#ffffff' },
        H:  { fill: '#f8fafc', hi: '#ffffff', rim: '#64748b', cov: 0.31, rad: 7.0,  txt: '#0f172a' }
      };
      const METALS = new Set(['Zr', 'Ti', 'Ga', 'Si']);

      function getStyle(sym) {
        return CPK_STYLE[sym] || { fill: '#64748b', hi: '#cbd5e1', rim: '#334155', cov: 1.0, rad: 11.0, txt: '#fff' };
      }

      function esc(s) {
        return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
      }

      function matMulVec(m, v) {
        return [
          m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
          m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
          m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
        ];
      }

      function rotateYawPitch(v, yaw, pitch) {
        const cy = Math.cos(yaw), sy = Math.sin(yaw);
        const cp = Math.cos(pitch), sp = Math.sin(pitch);
        const x1 = cy * v[0] + sy * v[2];
        const y1 = v[1];
        const z1 = -sy * v[0] + cy * v[2];
        return [
          x1,
          cp * y1 - sp * z1,
          sp * y1 + cp * z1
        ];
      }

      function buildSmoothSplinePath(pts) {
        if (pts.length < 2) return '';
        if (pts.length === 2) return 'M ' + pts[0].x.toFixed(1) + ',' + pts[0].y.toFixed(1) + ' L ' + pts[1].x.toFixed(1) + ',' + pts[1].y.toFixed(1);
        let yMin = Infinity, yMax = -Infinity;
        pts.forEach(p => { if (p.y < yMin) yMin = p.y; if (p.y > yMax) yMax = p.y; });
        yMin -= 6;
        yMax += 6;
        const out = [pts[0]];
        const segs = 12;
        for (let i = 0; i < pts.length - 1; i++) {
          const p0 = pts[Math.max(0, i - 1)];
          const p1 = pts[i];
          const p2 = pts[i + 1];
          const p3 = pts[Math.min(pts.length - 1, i + 2)];
          for (let s = 1; s <= segs; s++) {
            const t = s / segs, t2 = t * t, t3 = t2 * t;
            const x = 0.5 * ((2 * p1.x) + (-p0.x + p2.x) * t + (2 * p0.x - 5 * p1.x + 4 * p2.x - p3.x) * t2 + (-p0.x + 3 * p1.x - 3 * p2.x + p3.x) * t3);
            let y = 0.5 * ((2 * p1.y) + (-p0.y + p2.y) * t + (2 * p0.y - 5 * p1.y + 4 * p2.y - p3.y) * t2 + (-p0.y + 3 * p1.y - 3 * p2.y + p3.y) * t3);
            if (y < yMin) y = yMin;
            if (y > yMax) y = yMax;
            out.push({ x, y });
          }
        }
        return 'M ' + out.map(p => p.x.toFixed(1) + ',' + p.y.toFixed(1)).join(' L ');
      }

      function renderEnergyProfileSvg(rxn, activeStep) {
        const W = 470, H = 278;
        const x0 = 52, x1 = 446, y0 = 26, y1 = 232;
        const imgs = rxn.images || [];
        let maxE = Math.max(rxn.ea_lit || 0, rxn.ea_7n || 0, 20);
        let minE = 0;
        imgs.forEach(im => {
          if (im.energy_rel_kcal_mol > maxE) maxE = im.energy_rel_kcal_mol;
          if (im.energy_rel_kcal_mol < minE) minE = im.energy_rel_kcal_mol;
        });
        const yPadTop = Math.max(6, (maxE - minE) * 0.14);
        const yPadBot = Math.max(4, (maxE - minE) * 0.10);
        const loE = Math.floor((minE - yPadBot) / 5) * 5;
        const hiE = Math.ceil((maxE + yPadTop) / 5) * 5;

        const px = xi => x0 + xi * (x1 - x0);
        const py = e => y1 - ((e - loE) / Math.max(1e-6, hiE - loE)) * (y1 - y0);

        let grid = '';
        const nTicks = 5;
        for (let i = 0; i <= nTicks; i++) {
          const val = loE + (i / nTicks) * (hiE - loE);
          const yy = py(val);
          grid += '<line x1="' + x0 + '" y1="' + yy.toFixed(1) + '" x2="' + x1 + '" y2="' + yy.toFixed(1) + '" stroke="#e2e8f0" stroke-dasharray="3,3" />';
          grid += '<text x="' + (x0 - 7) + '" y="' + (yy + 3.5).toFixed(1) + '" text-anchor="end" font-size="9.5" fill="#475569">' + val.toFixed(0) + '</text>';
        }
        [0.0, 0.25, 0.5, 0.75, 1.0].forEach(t => {
          const xx = px(t);
          grid += '<line x1="' + xx.toFixed(1) + '" y1="' + y0 + '" x2="' + xx.toFixed(1) + '" y2="' + y1 + '" stroke="#f1f5f9" />';
          grid += '<text x="' + xx.toFixed(1) + '" y="' + (y1 + 14) + '" text-anchor="middle" font-size="9.2" fill="#475569">' + t.toFixed(2) + '</text>';
        });
        const yZero = py(0);
        grid += '<line x1="' + x0 + '" y1="' + yZero.toFixed(1) + '" x2="' + x1 + '" y2="' + yZero.toFixed(1) + '" stroke="#94a3b8" stroke-width="1.2" />';

        // Literature Ea reference line
        const yLit = py(rxn.ea_lit);
        let litLine = '<line x1="' + x0 + '" y1="' + yLit.toFixed(1) + '" x2="' + x1 + '" y2="' + yLit.toFixed(1) + '" stroke="#9333ea" stroke-width="1.4" stroke-dasharray="5,4" />';
        litLine += '<rect x="' + (x1 - 122) + '" y="' + (yLit - 15).toFixed(1) + '" width="120" height="14" rx="3" fill="#faf5ff" stroke="#d8b4fe" stroke-width="0.7" />';
        litLine += '<text x="' + (x1 - 62) + '" y="' + (yLit - 5).toFixed(1) + '" text-anchor="middle" font-size="8.4" font-weight="700" fill="#7e22ce">Lit Ea = ' + rxn.ea_lit.toFixed(2) + ' kcal/mol</text>';

        const pts = imgs.map(im => ({
          x: px(im.reaction_coord_norm),
          y: py(im.energy_rel_kcal_mol),
          im: im
        }));
        const splineD = buildSmoothSplinePath(pts);
        const areaD = splineD + ' L ' + pts[pts.length - 1].x.toFixed(1) + ',' + yZero.toFixed(1) + ' L ' + pts[0].x.toFixed(1) + ',' + yZero.toFixed(1) + ' Z';

        const actPt = pts[activeStep] || pts[0];
        const tsIdx = Math.min(pts.length - 1, Math.max(0, rxn.ts_image_index));
        const tsPt = pts[tsIdx];

        let nodesSvg = '';
        pts.forEach((p, idx) => {
          const isAct = (idx === activeStep);
          const isTs = (idx === tsIdx);
          let fill = '#ffffff';
          let stroke = '#0284c7';
          let rad = 4.5;
          if (isTs) {
            fill = '#7c3aed';
            stroke = '#5b21b6';
            rad = 6.2;
          }
          if (isAct) {
            stroke = '#d97706';
            if (!isTs) fill = '#f59e0b';
            rad = 7.0;
          }
          nodesSvg += '<g class="neb-plot-node" data-step="' + idx + '" style="cursor:pointer;">';
          nodesSvg += '<circle cx="' + p.x.toFixed(1) + '" cy="' + p.y.toFixed(1) + '" r="13" fill="transparent" />';
          if (isAct) {
            nodesSvg += '<circle cx="' + p.x.toFixed(1) + '" cy="' + p.y.toFixed(1) + '" r="10.5" fill="#fef3c7" opacity="0.85" />';
          }
          nodesSvg += '<circle cx="' + p.x.toFixed(1) + '" cy="' + p.y.toFixed(1) + '" r="' + rad + '" fill="' + fill + '" stroke="' + stroke + '" stroke-width="2" />';
          const lblY = p.y - 10 < y0 + 10 ? p.y + 16 : p.y - 9;
          let lblTxt = String(idx);
          if (idx === 0) lblTxt = 'IS(0)';
          else if (idx === pts.length - 1) lblTxt = 'FS(' + idx + ')';
          else if (isTs) lblTxt = 'TS(' + idx + ')';
          nodesSvg += '<text x="' + p.x.toFixed(1) + '" y="' + lblY.toFixed(1) + '" text-anchor="middle" font-size="8.4" font-weight="' + (isAct || isTs ? '800' : '600') + '" fill="' + (isAct ? '#b45309' : (isTs ? '#6d28d9' : '#334155')) + '">' + lblTxt + '</text>';
          nodesSvg += '</g>';
        });

        const crosshair = '<line x1="' + actPt.x.toFixed(1) + '" y1="' + y0 + '" x2="' + actPt.x.toFixed(1) + '" y2="' + y1 + '" stroke="#f59e0b" stroke-width="1.3" stroke-dasharray="3,2" />';
        const bannerTxt = 'Step ' + activeStep + '/' + (pts.length - 1) + '  |  xi=' + actPt.im.reaction_coord_norm.toFixed(3) + ' (s=' + actPt.im.reaction_coord_A.toFixed(2) + ' A)  |  dE=' + actPt.im.energy_rel_kcal_mol.toFixed(2) + ' kcal/mol';
        let tsTag = '';
        if (tsPt) {
          tsTag = '<text x="' + Math.min(x1 - 55, Math.max(x0 + 55, tsPt.x)).toFixed(1) + '" y="' + (tsPt.y - 20 < y0 + 8 ? tsPt.y + 26 : tsPt.y - 18).toFixed(1) + '" text-anchor="middle" font-size="8.8" font-weight="800" fill="#5b21b6">CI-TS: ' + tsPt.im.energy_rel_kcal_mol.toFixed(2) + ' kcal/mol</text>';
        }

        return '<svg viewBox="0 0 ' + W + ' ' + H + '" style="width:100%; height:auto; display:block; background:#ffffff; border:1px solid #cbd5e1; border-radius:6px;">' +
          '<defs><linearGradient id="nebAreaGrad" x1="0" y1="0" x2="0" y2="1"><stop offset="0%" stop-color="#0284c7" stop-opacity="0.22"/><stop offset="100%" stop-color="#0284c7" stop-opacity="0.02"/></linearGradient></defs>' +
          '<rect x="' + x0 + '" y="' + y0 + '" width="' + (x1 - x0) + '" height="' + (y1 - y0) + '" fill="#ffffff" />' +
          grid +
          '<path d="' + areaD + '" fill="url(#nebAreaGrad)" />' +
          '<path d="' + splineD + '" fill="none" stroke="#0284c7" stroke-width="2.5" stroke-linecap="round" />' +
          litLine + crosshair + tsTag + nodesSvg +
          '<rect x="' + (x0 + 4) + '" y="4" width="' + (x1 - x0 - 8) + '" height="17" rx="4" fill="#f1f5f9" stroke="#cbd5e1" stroke-width="0.8" />' +
          '<text x="' + ((x0 + x1) / 2).toFixed(1) + '" y="15.5" text-anchor="middle" font-size="9.0" font-weight="700" font-family="JetBrains Mono, monospace" fill="#0f172a">' + esc(bannerTxt) + '</text>' +
          '<text x="' + ((x0 + x1) / 2).toFixed(1) + '" y="266" text-anchor="middle" font-size="9.8" font-weight="700" fill="#334155">Normalized Reaction Coordinate &xi; = s / s_total (total s = ' + rxn.total_arc_length_A.toFixed(2) + ' A)</text>' +
          '<text x="15" y="' + ((y0 + y1) / 2).toFixed(1) + '" text-anchor="middle" font-size="9.8" font-weight="700" fill="#334155" transform="rotate(-90 15 ' + ((y0 + y1) / 2).toFixed(1) + ')">&Delta;E (kcal/mol)</text>' +
          '</svg>';
      }

      function renderStructure3dSvg(rxn, stepIdx, state) {
        const W = 460, H = 305;
        const img = (rxn.images && rxn.images[stepIdx]) || null;
        if (!img || !img.positions_A || img.positions_A.length === 0) {
          return '<svg viewBox="0 0 ' + W + ' ' + H + '" style="width:100%; height:auto; display:block; background:#0f172a; border-radius:6px;"><text x="230" y="152" text-anchor="middle" font-size="12" fill="#94a3b8">3D geometry snapshot unavailable</text></svg>';
        }
        const syms = rxn.symbols;
        const pos = img.positions_A;
        const basis = rxn.best_basis;
        const coreSet = new Set(rxn.reactive_indices || []);

        // Anchor center around midpoint/TS center across atoms so molecule stays centered across steps
        const tsImg = rxn.images[Math.min(rxn.images.length - 1, Math.max(0, rxn.ts_image_index))];
        let cx = 0, cy = 0, cz = 0;
        const tsPos = tsImg.positions_A;
        for (let i = 0; i < tsPos.length; i++) {
          cx += tsPos[i][0]; cy += tsPos[i][1]; cz += tsPos[i][2];
        }
        cx /= tsPos.length; cy /= tsPos.length; cz /= tsPos.length;

        const proj = pos.map((p, idx) => {
          const shifted = [p[0] - cx, p[1] - cy, p[2] - cz];
          const baseAligned = matMulVec(basis, shifted);
          const rot = rotateYawPitch(baseAligned, state.yaw, state.pitch);
          const camDist = Math.max(12.0, rxn.view_radius_A * 3.4);
          const persp = camDist / Math.max(3.5, camDist - rot[2]);
          const pxScale = (128.0 / Math.max(2.6, rxn.view_radius_A)) * state.zoom * persp;
          return {
            idx,
            sym: syms[idx],
            x3: p[0], y3: p[1], z3: p[2],
            rx: rot[0], ry: rot[1], rz: rot[2],
            sx: W * 0.5 + rot[0] * pxScale,
            sy: H * 0.5 - rot[1] * pxScale,
            persp,
            isCore: coreSet.has(idx)
          };
        });

        // Probe pair lookup
        const probeMap = new Map();
        const probes = rxn.bond_probe_pairs || [];
        probes.forEach(pr => {
          const a = Math.min(pr.i, pr.j), b = Math.max(pr.i, pr.j);
          probeMap.set(a + '_' + b, pr);
        });

        const drawItems = [];
        const nAtoms = pos.length;
        for (let i = 0; i < nAtoms; i++) {
          for (let j = i + 1; j < nAtoms; j++) {
            const dx = pos[i][0] - pos[j][0];
            const dy = pos[i][1] - pos[j][1];
            const dz = pos[i][2] - pos[j][2];
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
            const probe = probeMap.get(i + '_' + j);
            if (probe) {
              const isForm = probe.label.indexOf('forming') === 0;
              drawItems.push({
                kind: 'probe',
                z: 0.5 * (proj[i].rz + proj[j].rz) + 0.05,
                a: proj[i], b: proj[j],
                dist,
                isForm,
                label: probe.label
              });
              continue;
            }
            if (syms[i] === 'H' && syms[j] === 'H' && dist > 0.92) continue;
            const rSum = getStyle(syms[i]).cov + getStyle(syms[j]).cov;
            const hasMetal = METALS.has(syms[i]) || METALS.has(syms[j]);
            const tol = hasMetal ? 1.33 : 1.24;
            if (dist <= tol * rSum) {
              drawItems.push({
                kind: 'bond',
                z: 0.5 * (proj[i].rz + proj[j].rz),
                a: proj[i], b: proj[j],
                hasMetal
              });
            }
          }
        }

        proj.forEach(pt => {
          drawItems.push({ kind: 'atom', z: pt.rz, pt });
        });
        drawItems.sort((u, v) => u.z - v.z);

        let defs = '<defs>';
        Object.keys(CPK_STYLE).forEach(el => {
          const st = CPK_STYLE[el];
          defs += '<radialGradient id="cpk-' + state.uid + '-' + el + '" cx="34%" cy="30%" r="68%">' +
            '<stop offset="0%" stop-color="' + st.hi + '" />' +
            '<stop offset="55%" stop-color="' + st.fill + '" />' +
            '<stop offset="100%" stop-color="' + st.rim + '" />' +
            '</radialGradient>';
        });
        defs += '</defs>';

        let scene = '';
        const probeBadges = [];
        drawItems.forEach(it => {
          if (it.kind === 'bond') {
            const w = it.hasMetal ? 4.2 : 3.4;
            scene += '<line x1="' + it.a.sx.toFixed(1) + '" y1="' + it.a.sy.toFixed(1) + '" x2="' + it.b.sx.toFixed(1) + '" y2="' + it.b.sy.toFixed(1) + '" stroke="#475569" stroke-width="' + (w + 1.8) + '" stroke-linecap="round" />';
            scene += '<line x1="' + it.a.sx.toFixed(1) + '" y1="' + it.a.sy.toFixed(1) + '" x2="' + it.b.sx.toFixed(1) + '" y2="' + it.b.sy.toFixed(1) + '" stroke="#cbd5e1" stroke-width="' + w + '" stroke-linecap="round" opacity="0.88" />';
          } else if (it.kind === 'probe') {
            const col = it.isForm ? '#34d399' : '#fb923c';
            scene += '<line x1="' + it.a.sx.toFixed(1) + '" y1="' + it.a.sy.toFixed(1) + '" x2="' + it.b.sx.toFixed(1) + '" y2="' + it.b.sy.toFixed(1) + '" stroke="' + col + '" stroke-width="3.0" stroke-dasharray="5,4" stroke-linecap="round" opacity="0.95" />';
            probeBadges.push({
              mx: 0.5 * (it.a.sx + it.b.sx),
              my: 0.5 * (it.a.sy + it.b.sy),
              col,
              txt: it.dist.toFixed(2) + ' A'
            });
          } else {
            const p = it.pt;
            const st = getStyle(p.sym);
            const sr = st.rad * Math.min(1.3, Math.max(0.75, state.zoom * 0.92)) * p.persp;
            const gradId = CPK_STYLE[p.sym] ? ('url(#cpk-' + state.uid + '-' + p.sym + ')') : st.fill;
            if (p.isCore) {
              scene += '<circle cx="' + p.sx.toFixed(1) + '" cy="' + p.sy.toFixed(1) + '" r="' + (sr + 3.2).toFixed(1) + '" fill="none" stroke="#f59e0b" stroke-width="1.6" stroke-opacity="0.85" />';
            }
            scene += '<circle cx="' + p.sx.toFixed(1) + '" cy="' + p.sy.toFixed(1) + '" r="' + sr.toFixed(1) + '" fill="' + gradId + '" stroke="' + st.rim + '" stroke-width="0.9"><title>' + p.sym + p.idx + '</title></circle>';
            if (state.showLabels && (p.isCore || METALS.has(p.sym) || p.sym !== 'H')) {
              const showTxt = p.isCore ? (p.sym + p.idx) : p.sym;
              scene += '<text x="' + p.sx.toFixed(1) + '" y="' + (p.sy + 3.0).toFixed(1) + '" text-anchor="middle" font-size="' + (p.sym === 'H' ? '7.2' : '8.4') + '" font-weight="800" fill="' + st.txt + '" pointer-events="none">' + showTxt + '</text>';
            }
          }
        });

        probeBadges.forEach(bd => {
          scene += '<rect x="' + (bd.mx - 20).toFixed(1) + '" y="' + (bd.my - 8).toFixed(1) + '" width="40" height="14" rx="4" fill="#0f172a" fill-opacity="0.86" stroke="' + bd.col + '" stroke-width="1.0" pointer-events="none" />';
          scene += '<text x="' + bd.mx.toFixed(1) + '" y="' + (bd.my + 2).toFixed(1) + '" text-anchor="middle" font-size="8.0" font-weight="700" font-family="JetBrains Mono, monospace" fill="' + bd.col + '" pointer-events="none">' + bd.txt + '</text>';
        });

        const stepLabel = stepIdx === 0 ? 'Step 0: Initial State (IS)' : (stepIdx === rxn.ts_image_index ? 'Step ' + stepIdx + ': Climbing-Image TS' : (stepIdx === rxn.images.length - 1 ? 'Step ' + stepIdx + ': Final Adduct (FS)' : 'Step ' + stepIdx + ' / ' + (rxn.images.length - 1)));
        return '<svg viewBox="0 0 ' + W + ' ' + H + '" style="width:100%; height:auto; display:block; background:radial-gradient(circle at 50% 42%, #1e293b 0%, #0f172a 100%); border-radius:6px; border:1px solid #334155;">' +
          defs + scene +
          '<rect x="10" y="9" width="225" height="20" rx="4" fill="#0f172a" fill-opacity="0.82" stroke="#475569" stroke-width="0.8" />' +
          '<text x="18" y="22.5" font-size="9.5" font-weight="700" fill="#e2e8f0" pointer-events="none">' + esc(stepLabel) + ' | &Delta;E=' + img.energy_rel_kcal_mol.toFixed(2) + ' kcal/mol</text>' +
          '<rect x="10" y="' + (H - 24) + '" width="235" height="16" rx="3" fill="#0f172a" fill-opacity="0.72" />' +
          '<text x="16" y="' + (H - 12.5) + '" font-size="8.0" fill="#cbd5e1" pointer-events="none">Legend: <tspan fill="#fb923c" font-weight="700">-- Breaking</tspan>  <tspan fill="#34d399" font-weight="700">-- Forming</tspan>  <tspan fill="#fbbf24" font-weight="700">(O) Active Core</tspan></text>' +
          '</svg>';
      }

      function renderBondTableHtml(rxn, stepIdx) {
        const img = rxn.images[stepIdx] || {};
        const stepBonds = img.bond_lengths_A || {};
        const summaryBonds = rxn.bond_summary_A || {};
        const pairs = rxn.bond_probe_pairs || [];
        let rows = '';
        pairs.forEach(pr => {
          const sum = summaryBonds[pr.label] || {};
          const curVal = stepBonds[pr.label];
          const isForm = pr.label.indexOf('forming') === 0;
          const badge = isForm
            ? '<span style="color:#047857; font-weight:700;">Forming</span>'
            : '<span style="color:#c2410c; font-weight:700;">Breaking</span>';
          rows += '<tr>' +
            '<td style="padding:4px 6px;"><code>' + esc(pr.label) + '</code> (' + badge + ')</td>' +
            '<td style="padding:4px 6px; text-align:center; font-family:JetBrains Mono,monospace;">' + (sum.IS_A !== undefined ? sum.IS_A.toFixed(3) : '-') + '</td>' +
            '<td style="padding:4px 6px; text-align:center; font-family:JetBrains Mono,monospace; background:#fef3c7; font-weight:800; color:#92400e;">' + (curVal !== undefined ? curVal.toFixed(3) : '-') + '</td>' +
            '<td style="padding:4px 6px; text-align:center; font-family:JetBrains Mono,monospace;">' + (sum.TS_A !== undefined ? sum.TS_A.toFixed(3) : '-') + '</td>' +
            '<td style="padding:4px 6px; text-align:center; font-family:JetBrains Mono,monospace;">' + (sum.FS_A !== undefined ? sum.FS_A.toFixed(3) : '-') + '</td>' +
            '</tr>';
        });
        return '<div style="font-size:0.77rem; margin-bottom:6px; color:#334155;"><strong>반응식:</strong> <code>' + esc(rxn.equation) + '</code> &bull; <strong>기구:</strong> ' + esc(rxn.mechanism) + '</div>' +
          '<div style="overflow-x:auto;"><table style="width:100%; border-collapse:collapse; font-size:0.74rem; background:#fff; border:1px solid #e2e8f0; border-radius:4px;">' +
          '<thead><tr style="background:#f1f5f9; color:#334155;"><th style="padding:4px 6px; text-align:left;">활성 결합 (A)</th><th style="padding:4px 6px;">IS(0)</th><th style="padding:4px 6px; background:#fde68a;">현재(Step ' + stepIdx + ')</th><th style="padding:4px 6px;">TS(' + rxn.ts_image_index + ')</th><th style="padding:4px 6px;">FS(' + (rxn.images.length - 1) + ')</th></tr></thead>' +
          '<tbody>' + rows + '</tbody></table></div>';
      }

      function mountWidget(el, reactions) {
        const eid = el.getAttribute('data-explorer-id');
        const byId = {};
        reactions.forEach(r => { byId[r.reaction_id] = r; });
        const initId = byId[el.getAttribute('data-default-rxn')] ? el.getAttribute('data-default-rxn') : reactions[0].reaction_id;

        const state = {
          uid: eid,
          rxnId: initId,
          step: byId[initId].ts_image_index || 0,
          yaw: 0.0,
          pitch: 0.0,
          zoom: 1.0,
          showLabels: true,
          timer: null
        };

        const tabsHost = document.getElementById('neb-tabs-' + eid);
        const energyHost = document.getElementById('neb-energy-svg-host-' + eid);
        const bondHost = document.getElementById('neb-bond-table-host-' + eid);
        const structHost = document.getElementById('neb-struct-svg-host-' + eid);
        const stripHost = document.getElementById('neb-step-strip-' + eid);
        const slider = document.getElementById('neb-slider-' + eid);
        const playBtn = document.getElementById('neb-play-btn-' + eid);
        const badgeEl = document.getElementById('neb-badge-' + eid);
        const curveMetaEl = document.getElementById('neb-curve-meta-' + eid);

        function stopPlay() {
          if (state.timer) {
            clearInterval(state.timer);
            state.timer = null;
          }
          if (playBtn) playBtn.textContent = '자동 재생';
        }

        function draw() {
          const rxn = byId[state.rxnId];
          const maxIdx = rxn.images.length - 1;
          if (state.step < 0) state.step = 0;
          if (state.step > maxIdx) state.step = maxIdx;
          if (slider) {
            slider.max = String(maxIdx);
            slider.value = String(state.step);
          }

          if (tabsHost) {
            tabsHost.innerHTML = reactions.map(r => {
              const active = (r.reaction_id === state.rxnId);
              const bg = active ? '#4f46e5' : '#f1f5f9';
              const fg = active ? '#ffffff' : '#334155';
              const br = active ? '#4338ca' : '#cbd5e1';
              return '<button type="button" data-rxn="' + r.reaction_id + '" style="padding:5px 10px; font-size:0.77rem; font-weight:700; border-radius:999px; border:1px solid ' + br + '; background:' + bg + '; color:' + fg + '; cursor:pointer;">' + esc(r.short) + ' (' + r.ea_7n.toFixed(1) + ')</button>';
            }).join('');
            tabsHost.querySelectorAll('button[data-rxn]').forEach(btn => {
              btn.addEventListener('click', () => {
                stopPlay();
                state.rxnId = btn.getAttribute('data-rxn');
                const nr = byId[state.rxnId];
                state.step = nr.ts_image_index || 0;
                state.yaw = 0.0;
                state.pitch = 0.0;
                state.zoom = 1.0;
                draw();
              });
            });
          }

          if (badgeEl) {
            badgeEl.textContent = rxn.short + ' | 7net Ea=' + rxn.ea_7n.toFixed(2) + ' vs Lit ' + rxn.ea_lit.toFixed(2) + ' kcal/mol | n_imag=' + rxn.n_imag;
          }
          if (curveMetaEl) {
            curveMetaEl.textContent = 'TS=Img ' + rxn.ts_image_index + ' | fmax=' + rxn.fmax.toFixed(3) + ' eV/A';
          }

          if (energyHost) {
            energyHost.innerHTML = renderEnergyProfileSvg(rxn, state.step);
            energyHost.querySelectorAll('.neb-plot-node').forEach(node => {
              node.addEventListener('click', () => {
                stopPlay();
                state.step = parseInt(node.getAttribute('data-step'), 10);
                draw();
              });
            });
          }

          if (bondHost) bondHost.innerHTML = renderBondTableHtml(rxn, state.step);
          if (structHost) structHost.innerHTML = renderStructure3dSvg(rxn, state.step, state);

          if (stripHost) {
            stripHost.innerHTML = rxn.images.map((im, idx) => {
              const isAct = (idx === state.step);
              const isTs = (idx === rxn.ts_image_index);
              let tag = 'Step ' + idx;
              if (idx === 0) tag = '0 (IS)';
              else if (isTs) tag = idx + ' (TS)';
              else if (idx === maxIdx) tag = idx + ' (FS)';
              const br = isAct ? '#d97706' : (isTs ? '#7c3aed' : '#cbd5e1');
              const bg = isAct ? '#fef3c7' : (isTs ? '#f5f3ff' : '#ffffff');
              return '<button type="button" data-step="' + idx + '" style="padding:4px 2px; border-radius:5px; border:' + (isAct ? '2px' : '1px') + ' solid ' + br + '; background:' + bg + '; cursor:pointer; text-align:center; line-height:1.25;">' +
                '<div style="font-size:0.7rem; font-weight:800; color:#0f172a;">' + tag + '</div>' +
                '<div style="font-size:0.65rem; font-family:JetBrains Mono,monospace; color:#475569;">&xi;=' + im.reaction_coord_norm.toFixed(2) + '</div>' +
                '<div style="font-size:0.67rem; font-weight:700; color:' + (isTs ? '#6d28d9' : '#0369a1') + ';">' + im.energy_rel_kcal_mol.toFixed(1) + '</div>' +
                '</button>';
            }).join('');
            stripHost.querySelectorAll('button[data-step]').forEach(btn => {
              btn.addEventListener('click', () => {
                stopPlay();
                state.step = parseInt(btn.getAttribute('data-step'), 10);
                draw();
              });
            });
          }
        }

        if (slider) {
          slider.addEventListener('input', () => {
            stopPlay();
            state.step = parseInt(slider.value, 10);
            draw();
          });
        }

        el.querySelectorAll('.neb-ctrl-btn').forEach(btn => {
          btn.addEventListener('click', () => {
            const act = btn.getAttribute('data-act');
            const maxIdx = byId[state.rxnId].images.length - 1;
            if (act === 'prev') {
              stopPlay();
              state.step = (state.step - 1 + maxIdx + 1) % (maxIdx + 1);
              draw();
            } else if (act === 'next') {
              stopPlay();
              state.step = (state.step + 1) % (maxIdx + 1);
              draw();
            } else if (act === 'jump-ts') {
              stopPlay();
              state.step = byId[state.rxnId].ts_image_index || 0;
              draw();
            } else if (act === 'zoom-in') {
              state.zoom = Math.min(2.2, state.zoom * 1.15);
              draw();
            } else if (act === 'zoom-out') {
              state.zoom = Math.max(0.55, state.zoom / 1.15);
              draw();
            } else if (act === 'reset-view') {
              state.yaw = 0.0;
              state.pitch = 0.0;
              state.zoom = 1.0;
              draw();
            } else if (act === 'toggle-labels') {
              state.showLabels = !state.showLabels;
              draw();
            } else if (act === 'play') {
              if (state.timer) {
                stopPlay();
              } else {
                if (playBtn) playBtn.textContent = '일시 정지';
                state.timer = setInterval(() => {
                  const mx = byId[state.rxnId].images.length - 1;
                  state.step = (state.step + 1) % (mx + 1);
                  draw();
                }, 900);
              }
            }
          });
        });

        // 3D pointer orbit rotation on structHost
        if (structHost) {
          let dragging = false, lastX = 0, lastY = 0;
          structHost.addEventListener('pointerdown', ev => {
            dragging = true;
            lastX = ev.clientX;
            lastY = ev.clientY;
            structHost.style.cursor = 'grabbing';
            if (structHost.setPointerCapture) structHost.setPointerCapture(ev.pointerId);
          });
          structHost.addEventListener('pointermove', ev => {
            if (!dragging) return;
            const dx = ev.clientX - lastX;
            const dy = ev.clientY - lastY;
            lastX = ev.clientX;
            lastY = ev.clientY;
            state.yaw += dx * 0.012;
            state.pitch += dy * 0.012;
            structHost.innerHTML = renderStructure3dSvg(byId[state.rxnId], state.step, state);
          });
          const endDrag = () => {
            dragging = false;
            structHost.style.cursor = 'grab';
          };
          structHost.addEventListener('pointerup', endDrag);
          structHost.addEventListener('pointercancel', endDrag);
        }

        draw();
      }

      document.addEventListener('DOMContentLoaded', function() {
        const scriptEl = document.getElementById('neb-trajectory-payload');
        if (!scriptEl) return;
        let payload = null;
        try {
          payload = JSON.parse(scriptEl.textContent || '{}');
        } catch (e) {
          console.warn('NEB payload parse warning:', e);
          return;
        }
        const reactions = payload.reactions || [];
        if (reactions.length === 0) return;
        document.querySelectorAll('.neb-explorer-widget').forEach(el => mountWidget(el, reactions));
      });
    })();
