"""Introduction Tab Builder (``소개``) for SAPIS.

Generates the standalone HTML content for the model introduction tab
(``#tab-content-intro``, button label ``소개``), documenting the official
model name ``SAPIS``, its full name, Latin etymology, and 2-Track scoring
architecture.
"""

from __future__ import annotations


def build_intro_html() -> str:
    """Build the ``소개`` tab HTML fragment for SAPIS identity, etymology, and architecture."""
    return """
      <section class="card" style="background:linear-gradient(135deg, #f8fafc 0%, #eff6ff 100%); border:1px solid #bfdbfe; margin-bottom:18px;">
        <div style="display:flex; flex-wrap:wrap; justify-content:space-between; align-items:flex-start; gap:16px;">
          <div style="flex:1 1 540px;">
            <div style="display:inline-block; font-size:0.74rem; font-weight:700; letter-spacing:0.06em; text-transform:uppercase; color:#1d4ed8; background:#dbeafe; border:1px solid #93c5fd; padding:3px 10px; border-radius:999px; margin-bottom:10px;">
              Unified Semiconductor SA Model
            </div>
            <h2 style="margin:0 0 8px 0; font-size:1.75rem; font-weight:800; color:#0f172a; letter-spacing:-0.02em;">
              SAPIS
            </h2>
            <p style="margin:0 0 10px 0; font-size:1.02rem; font-weight:700; color:#0369a1; font-family:'JetBrains Mono', monospace;">
              Synthetic Accessibility of Precursors &amp; Inhibitors for Semiconductors
            </p>
            <p style="margin:0; font-size:0.9rem; color:#334155; line-height:1.7;">
              <strong>SAPIS</strong>는 반도체 박막 증착용 <strong>ALD/CVD 프리커서</strong>(전이금속 및 주족·준금속 포함)와
              선택적 원자층 증착용 <strong>AS-ALD 인히비터</strong>의 합성 가능성을 단일 기준으로 통합 평가하기 위해 설계된
              <strong>2-Track 물리·화학 정보학 하이브리드 모델</strong>입니다.
            </p>
          </div>
          <div style="flex:0 0 auto; min-width:240px; background:#ffffff; border:1px solid #cbd5e1; border-radius:10px; padding:14px 16px;">
            <div style="font-size:0.75rem; font-weight:700; color:#64748b; margin-bottom:6px;">핵심 정체성 요약</div>
            <ul style="margin:0; padding-left:18px; font-size:0.84rem; color:#1e293b; line-height:1.7;">
              <li><strong>공식 명칭:</strong> SAPIS</li>
              <li><strong>평가 대상:</strong> 반도체 프리커서 &amp; 인히비터</li>
              <li><strong>점수 범위:</strong> 1.0 (용이) ~ 10.0 (난해)</li>
              <li><strong>합성 가능 기준:</strong> $\\le 3.5$ (Pass)</li>
            </ul>
          </div>
        </div>
      </section>

      <section class="card" style="margin-bottom:18px;">
        <h3 style="margin:0 0 12px 0; font-size:1.12rem; font-weight:800; color:#0f172a;">
          1. 철자 구성 (Acronym Breakdown)
        </h3>
        <p style="margin:0 0 14px 0; font-size:0.88rem; color:#475569; line-height:1.65;">
          <strong>SAPIS</strong>의 다섯 글자는 불필요한 연결어나 보조 알파벳 없이 본 모델이 다루는 4대 핵심 요소만으로 구성됩니다.
        </p>
        <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(210px, 1fr)); gap:12px;">
          <div style="background:#eff6ff; border:1px solid #bfdbfe; border-radius:10px; padding:14px;">
            <div style="font-size:1.25rem; font-weight:800; color:#1d4ed8; font-family:'JetBrains Mono', monospace;">SA</div>
            <div style="font-size:0.88rem; font-weight:700; color:#0f172a; margin:4px 0;">Synthetic Accessibility</div>
            <div style="font-size:0.81rem; color:#475569; line-height:1.6;">
              분자 조각 빈도와 리간드 분리(<code>ligand_sa</code>)를 바탕으로 한 정량적 합성 난이도 평가 체계
            </div>
          </div>
          <div style="background:#f0fdfa; border:1px solid #99f6e4; border-radius:10px; padding:14px;">
            <div style="font-size:1.25rem; font-weight:800; color:#0f766e; font-family:'JetBrains Mono', monospace;">P</div>
            <div style="font-size:0.88rem; font-weight:700; color:#0f172a; margin:4px 0;">Precursors</div>
            <div style="font-size:0.81rem; color:#475569; line-height:1.6;">
              Hf, Zr, Ti 등 전이금속과 Si, Ge, Al, B 등 주족·준금속을 아우르는 반도체 ALD/CVD 프리커서
            </div>
          </div>
          <div style="background:#fefce8; border:1px solid #fde68a; border-radius:10px; padding:14px;">
            <div style="font-size:1.25rem; font-weight:800; color:#b45309; font-family:'JetBrains Mono', monospace;">I</div>
            <div style="font-size:0.88rem; font-weight:700; color:#0f172a; margin:4px 0;">Inhibitors</div>
            <div style="font-size:0.81rem; color:#475569; line-height:1.6;">
              영역 선택적 원자층 증착(AS-ALD)에서 특정 표면을 패시베이션하는 유기 표면 차단제
            </div>
          </div>
          <div style="background:#f8fafc; border:1px solid #cbd5e1; border-radius:10px; padding:14px;">
            <div style="font-size:1.25rem; font-weight:800; color:#334155; font-family:'JetBrains Mono', monospace;">S</div>
            <div style="font-size:0.88rem; font-weight:700; color:#0f172a; margin:4px 0;">Semiconductors</div>
            <div style="font-size:0.81rem; color:#475569; line-height:1.6;">
              반도체 공정 특화 도메인(고순도·열안정성·휘발성 및 3D 배위 물리 법칙 검증 포함)
            </div>
          </div>
        </div>
      </section>

      <section class="card" style="margin-bottom:18px;">
        <h3 style="margin:0 0 10px 0; font-size:1.12rem; font-weight:800; color:#0f172a;">
          2. 어원 및 명명 철학 (Etymology &amp; Rationale)
        </h3>
        <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(290px, 1fr)); gap:14px; align-items:stretch;">
          <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px; padding:16px;">
            <div style="font-size:0.78rem; font-weight:700; color:#0369a1; margin-bottom:6px;">라틴어 어원</div>
            <p style="margin:0 0 8px 0; font-size:0.9rem; color:#0f172a; line-height:1.7;">
              <strong>SAPIS</strong>는 <em>Homo sapiens</em>의 어근인 라틴어 동사
              <strong><em>sapere</em></strong> (알다, 맛을 보다, 통찰하다, 분별하다)의
              2인칭 단수 직설법 현재형인 <strong><em>sapis</em></strong>(&ldquo;너는 분별한다&rdquo;)와 철자가 일치합니다.
            </p>
            <p style="margin:0; font-size:0.84rem; color:#475569; line-height:1.65;">
              단순히 데이터셋 빈도표만 기계적으로 조회하는 수준을 넘어, 화학적 합성 경로와 3D 물리적 안정성을 함께 읽어내어
              공정 후보 물질의 옥석을 <strong>분별한다</strong>는 모델의 지향점을 담았습니다.
            </p>
          </div>
          <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px; padding:16px;">
            <div style="font-size:0.78rem; font-weight:700; color:#0f766e; margin-bottom:6px;">설계 일관성</div>
            <p style="margin:0 0 8px 0; font-size:0.9rem; color:#0f172a; line-height:1.7;">
              기존 2D 유기 SA 모델들은 <code>TiEt4</code>나 <code>ZrMe4</code>처럼 열역학적으로 불안정한
              금속 알킬 착체도 조각 빈도만 보고 합성 가능(Pass)으로 오탐하는 한계가 있었습니다.
            </p>
            <p style="margin:0; font-size:0.84rem; color:#475569; line-height:1.65;">
              <strong>SAPIS</strong>는 Track 1의 리간드 합성성(<code>ligand_sa</code>)과
              Track 2의 SevenNet-0 기반 3D 물리 게이트($\\mathrm{BDE}_{\\min}$, $\\%V_{\\mathrm{bur}}$, $E_{\\mathrm{coh}}$, Hydride)를 결합해
              실제 합성·공정 가능성을 물리적으로 분별합니다.
            </p>
          </div>
        </div>
      </section>

      <section class="card">
        <h3 style="margin:0 0 10px 0; font-size:1.12rem; font-weight:800; color:#0f172a;">
          3. SAPIS 통합 스코어링 구조 요약
        </h3>
        <div style="background:#f8fafc; border:1px solid #cbd5e1; border-radius:8px; padding:12px 16px; margin-bottom:12px; font-size:0.9rem; overflow-x:auto;">
          $$\\mathrm{SAPIS} = \\operatorname{clip}\\!\\Big(\\text{Track 1}\\;(\\texttt{ligand\\_sa}) + \\text{Track 2}\\;(\\Delta_{\\text{labile}} + \\Delta_{\\text{crowd}} + \\Delta_{\\text{olig}} + \\Delta_{\\text{hydride}}),\\; 1.0,\\; 10.0\\Big)$$
        </div>
        <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(260px, 1fr)); gap:12px; font-size:0.84rem; color:#334155;">
          <div style="border-left:3px solid #2563eb; padding-left:12px;">
            <strong style="color:#0f172a;">Track 1 — 화학적 합성 난이도 (<code>ligand_sa</code>)</strong>
            <p style="margin:4px 0 0 0; line-height:1.65;">
              프리커서는 금속/준금속 중심과 프로리간드($\\mathrm{HL}_i$)를 분리한 뒤 각 리간드의
              $\\min(\\text{BR-SAScore},\\, \\text{SAScore})$ 최댓값으로 평가하고, 유기 인히비터는 분자 전체에 대해 동일 기준을 적용합니다.
              이때 <code>BR-SAScore</code>의 미등재 조각 패널티(<code>-6.0</code>)로 인해 상용 시약(예: <code>Hacac</code>, <code>5.20</code>)이 과대 평가되는 편향을,
              동일한 <code>1.0~10.0</code> 스케일 및 조각 빈도 체계를 공유하는 <code>SAScore</code>(<code>2.17</code>)와의 최솟값 결합으로 보완합니다.
            </p>
          </div>
          <div style="border-left:3px solid #0d9488; padding-left:12px;">
            <strong style="color:#0f172a;">Track 2 — 3D 물리 안정성 페널티 (SevenNet-0)</strong>
            <p style="margin:4px 0 0 0; line-height:1.65;">
              최소 결합 해리 에너지($\\mathrm{BDE}_{\\min}$), 입체 차폐율($\\%V_{\\mathrm{bur}}$),
              응집 에너지($E_{\\mathrm{coh}}$), 불안정 하이드라이드($\\Delta_{\\text{hydride}}$) 게이트를 합산해
              합성·저장 불가능한 불안정 후보를 걸러냅니다.
            </p>
          </div>
        </div>
      </section>
    """
