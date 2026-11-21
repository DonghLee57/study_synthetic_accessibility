    const KATEX_RENDER_OPTIONS = {
      delimiters: [
        {left: '$$', right: '$$', display: true},
        {left: '$', right: '$', display: false},
        {left: '\\[', right: '\\]', display: true},
        {left: '\\(', right: '\\)', display: false}
      ],
      throwOnError: false
    };
    const renderedKatexTabs = new Set();

    function renderKatexInContainer(tabKey, targetEl) {
      if (typeof renderMathInElement !== 'function' || !targetEl) return;
      if (tabKey && renderedKatexTabs.has(tabKey)) return;
      renderMathInElement(targetEl, KATEX_RENDER_OPTIONS);
      if (tabKey) renderedKatexTabs.add(tabKey);
    }

    document.addEventListener("DOMContentLoaded", function() {
      // Render math only in visible header & initial intro tab (defer hidden tabs)
      renderKatexInContainer('header', document.querySelector('header'));
      renderKatexInContainer('intro', document.getElementById('tab-content-intro'));
      initScrollSpy();

      // Check URL hash on direct link entry
      if (window.location.hash) {
        const targetId = window.location.hash.substring(1);
        const targetEl = document.getElementById(targetId);
        const contentTechnical = document.getElementById('tab-content-technical');
        if (targetEl && contentTechnical && contentTechnical.contains(targetEl)) {
          switchViewTab('technical', { skipScrollTop: true });
          requestAnimationFrame(() => {
            targetEl.scrollIntoView({ behavior: 'auto' });
          });
        }
      }
    });

    // Individual Table Toggle
    function togglePanel(panelId, btnId, closedText, openText) {
      const panel = document.getElementById(panelId);
      const btn = document.getElementById(btnId);
      if (!panel || !btn) return;

      const isOpen = panel.classList.contains('open');
      if (isOpen) {
        panel.classList.remove('open');
        btn.classList.remove('is-open');
        btn.textContent = closedText;
      } else {
        panel.classList.add('open');
        btn.classList.add('is-open');
        btn.textContent = openText;
      }
      updateGlobalToggleBtnState();
    }

    // Global Expand/Collapse All Tables from Sidebar
    const tableConfigs = [
      { panel: 'panel-lit', btn: 'btn-panel-lit', closed: '문헌 검증 전체 테이블 보기 (25건)', open: '문헌 검증 테이블 닫기 (숨김)' },
      { panel: 'panel-matrix', btn: 'btn-panel-matrix', closed: '16개 모델 비교 테이블 보기', open: '비교 테이블 닫기 (숨김)' },
      
      { panel: 'panel-bottleneck', btn: 'btn-panel-bottleneck', closed: '3대 병목 및 해결책 테이블 보기', open: '병목 비교 테이블 닫기 (숨김)' },
      { panel: 'panel-vasp', btn: 'btn-panel-vasp', closed: 'VASP 자동화 프레임워크 테이블 보기 (10건)', open: 'VASP 프레임워크 테이블 닫기 (숨김)' },
      { panel: 'panel-precdb', btn: 'btn-panel-precdb', closed: '상용 프리커서 DB 테이블 보기 (27종)', open: '상용 프리커서 DB 닫기 (숨김)' },
      { panel: 'panel-bbdb', btn: 'btn-panel-bbdb', closed: '출발물질·시약 DB 테이블 보기 (29종)', open: '출발물질·시약 DB 닫기 (숨김)' },
      { panel: 'panel-rxndb', btn: 'btn-panel-rxndb', closed: '공인 합성 반응 DB 테이블 보기 (12종)', open: '합성 반응 DB 닫기 (숨김)' },
      { panel: 'panel-expdb', btn: 'btn-panel-expdb', closed: '확장 프리커서 DB 테이블 보기 (205종)', open: '확장 프리커서 DB 닫기 (숨김)' },
      { panel: 'panel-inhibdb', btn: 'btn-panel-inhibdb', closed: 'AS-ALD 억제제 DB 테이블 보기 (18종)', open: 'AS-ALD 억제제 DB 닫기 (숨김)' },
      { panel: 'panel-inversedb', btn: 'btn-panel-inversedb', closed: '물성 지표 & 연산 매핑 테이블 보기 (7종)', open: '매핑 테이블 닫기 (숨김)' }
    ];

    function toggleAllTables() {
      const anyClosed = tableConfigs.some(cfg => !document.getElementById(cfg.panel).classList.contains('open'));
      tableConfigs.forEach(cfg => {
        const panel = document.getElementById(cfg.panel);
        const btn = document.getElementById(cfg.btn);
        if (!panel || !btn) return;
        if (anyClosed) {
          panel.classList.add('open');
          btn.classList.add('is-open');
          btn.textContent = cfg.open;
        } else {
          panel.classList.remove('open');
          btn.classList.remove('is-open');
          btn.textContent = cfg.closed;
        }
      });
      updateGlobalToggleBtnState();
    }

    function updateGlobalToggleBtnState() {
      const globalBtn = document.getElementById('btnToggleAllTables');
      const anyClosed = tableConfigs.some(cfg => !document.getElementById(cfg.panel).classList.contains('open'));
      globalBtn.textContent = anyClosed ? '전체 상세 테이블 일괄 펼치기' : '전체 상세 테이블 일괄 접기';
    }

    // Filter rows inside Section 1 Table
    function filterPapers(gen, event) {
      const buttons = document.querySelectorAll('.filter-btn');
      buttons.forEach(btn => btn.classList.remove('active'));
      if (event && event.target) {
        event.target.classList.add('active');
      }

      const rows = document.querySelectorAll('.lit-row');
      rows.forEach(row => {
        const rowGen = row.getAttribute('data-gen');
        const downloaded = row.getAttribute('data-downloaded') === 'true';

        if (gen === 'all') {
          row.style.display = '';
        } else if (gen === 'downloaded') {
          row.style.display = downloaded ? '' : 'none';
        } else {
          row.style.display = rowGen === gen ? '' : 'none';
        }
      });
    }

    function searchPapers() {
      const query = document.getElementById('paperSearch').value.toLowerCase();
      const rows = document.querySelectorAll('.lit-row');
      rows.forEach(row => {
        const text = row.getAttribute('data-search');
        row.style.display = text.includes(query) ? '' : 'none';
      });
    }

    // Filter & Search for Section 11 Precursor DB
    function filterPrecursors(cat, event) {
      const panel = document.getElementById('panel-precdb');
      const buttons = panel.querySelectorAll('.filter-btn');
      buttons.forEach(btn => btn.classList.remove('active'));
      if (event && event.target) {
        event.target.classList.add('active');
      }

      const rows = document.querySelectorAll('.prec-row');
      rows.forEach(row => {
        const rowCat = row.getAttribute('data-cat');
        if (cat === 'all' || rowCat === cat) {
          row.style.display = '';
        } else {
          row.style.display = 'none';
        }
      });
    }

    function searchPrecursors() {
      const raw = document.getElementById('precSearch').value.trim();
      const query = raw.toLowerCase();
      const metalSym = typeof resolveMetalQuery === 'function' ? resolveMetalQuery(raw) : null;
      const rows = document.querySelectorAll('.prec-row');
      rows.forEach(row => {
        const text = row.getAttribute('data-search') || '';
        const cellElem = row.children[3] ? row.children[3].textContent.trim() : '';
        const matched = !query ? true : (metalSym ? cellElem === metalSym : text.includes(query));
        row.style.display = matched ? '' : 'none';
      });
    }

    function searchExpandedPrecursors() {
      const input = document.getElementById('expPrecSearch');
      if (!input) return;
      const raw = input.value.trim();
      const query = raw.toLowerCase();
      const metalSym = typeof resolveMetalQuery === 'function' ? resolveMetalQuery(raw) : null;
      const rows = document.querySelectorAll('.exp-prec-row');
      rows.forEach(row => {
        const text = row.getAttribute('data-search') || '';
        const cellElem = row.children[3] ? row.children[3].textContent.trim() : '';
        const matched = !query ? true : (metalSym ? cellElem === metalSym : text.includes(query));
        row.style.display = matched ? '' : 'none';
      });
    }

    function searchInhibitors() {
      const input = document.getElementById('inhibSearch');
      if (!input) return;
      const query = input.value.toLowerCase();
      const rows = document.querySelectorAll('.inhib-row');
      rows.forEach(row => {
        const text = row.getAttribute('data-search') || '';
        row.style.display = text.includes(query) ? '' : 'none';
      });
    }

    function filterBenchmark(type, event) {
      const panel = document.getElementById('panel-benchmark');
      if (panel) {
        const buttons = panel.querySelectorAll('.filter-btn');
        buttons.forEach(btn => btn.classList.remove('active'));
      }
      if (event && event.target) {
        event.target.classList.add('active');
      }
      const rows = document.querySelectorAll('.benchmark-row');
      rows.forEach(row => {
        const rowType = row.getAttribute('data-type');
        if (type === 'all' || rowType === type) {
          row.style.display = '';
        } else {
          row.style.display = 'none';
        }
      });
    }

    function searchBenchmark() {
      const input = document.getElementById('benchmarkSearch');
      if (!input) return;
      const query = input.value.toLowerCase();
      const rows = document.querySelectorAll('.benchmark-row');
      rows.forEach(row => {
        const text = row.getAttribute('data-search') || '';
        row.style.display = text.includes(query) ? '' : 'none';
      });
    }

    // Diagram Popup Modal: Direct Click to Enlarge
    let currentDiagramZoom = 1.0;



    function updateDiagramTransform() {

      const container = document.getElementById('modalSvgContainer');

      const badge = document.getElementById('modalZoomLevel');

      if (!container) return;

      container.style.transform = `scale(${currentDiagramZoom})`;

      if (badge) {

        badge.textContent = `${Math.round(currentDiagramZoom * 100)}%`;

      }

    }



    function zoomDiagram(delta) {

      currentDiagramZoom = Math.max(0.3, Math.min(3.5, currentDiagramZoom + delta));

      updateDiagramTransform();

    }



    function resetDiagramZoom() {

      currentDiagramZoom = 1.0;

      updateDiagramTransform();

    }



    function fitDiagramToWindow() {

      const body = document.getElementById('modalDiagramBody');

      const container = document.getElementById('modalSvgContainer');

      const svg = container ? container.querySelector('svg') : null;

      if (!body || !svg) return;



      // Extract native viewBox dimensions

      const viewBox = svg.viewBox && svg.viewBox.baseVal ? svg.viewBox.baseVal : null;

      let svgW = viewBox && viewBox.width > 0 ? viewBox.width : svg.getBoundingClientRect().width;

      let svgH = viewBox && viewBox.height > 0 ? viewBox.height : svg.getBoundingClientRect().height;



      const availW = body.clientWidth - 48;

      const availH = body.clientHeight - 48;



      if (svgW > 0 && svgH > 0 && availW > 0 && availH > 0) {

        const scaleX = availW / svgW;

        const scaleY = availH / svgH;

        // Fit within viewport with 92% fill factor

        currentDiagramZoom = Math.min(2.5, Math.max(0.4, Math.min(scaleX, scaleY) * 0.95));

      } else {

        currentDiagramZoom = 1.0;

      }

      updateDiagramTransform();

    }



    async function openDiagramModal(sourceContainerId, title) {
      const modal = document.getElementById('diagramModal');
      const modalTitle = document.getElementById('modalDiagramTitle');
      const modalContainer = document.getElementById('modalSvgContainer');
      const sourceEl = document.getElementById(sourceContainerId);
      if (!modal || !sourceEl) return;

      // Guarantee lazy Mermaid diagram is rendered before cloning into modal
      if (sourceEl.classList.contains('mermaid') && sourceEl.getAttribute('data-processed') !== 'true') {
        await renderSingleMermaidElement(sourceEl, 'modal');
      } else {
        const pendingChild = sourceEl.querySelector('.mermaid:not([data-processed="true"])');
        if (pendingChild) {
          await renderSingleMermaidElement(pendingChild, 'modal');
        }
      }

      modalTitle.textContent = title || '다이어그램 확대 팝업';
      modalContainer.innerHTML = '';
      
      const sourceSvg = sourceEl.querySelector('svg');
      if (sourceSvg) {
        const clonedSvg = sourceSvg.cloneNode(true);
        modalContainer.appendChild(clonedSvg);
        const viewBox = clonedSvg.viewBox && clonedSvg.viewBox.baseVal ? clonedSvg.viewBox.baseVal : null;
        if (viewBox && viewBox.width > 0) {
          clonedSvg.style.width = Math.max(viewBox.width, 1050) + 'px';
          clonedSvg.style.height = 'auto';
          clonedSvg.style.maxWidth = 'none';
        } else {
          clonedSvg.style.width = '100%';
          clonedSvg.style.height = 'auto';
        }
      } else {
        modalContainer.innerHTML = sourceEl.innerHTML;
      }

      currentDiagramZoom = 1.0;
      updateDiagramTransform();

      modal.classList.add('active');
      document.body.style.overflow = 'hidden';

      // Automatically scale up to comfortably fit the modal window
      setTimeout(() => {
        fitDiagramToWindow();
      }, 40);
    }



    function closeDiagramModal() {

      const modal = document.getElementById('diagramModal');

      if (modal) {

        modal.classList.remove('active');

        document.body.style.overflow = '';

      }

    }

    function closeDiagramModalOnOverlay(event) {
      if (event.target && event.target.id === 'diagramModal') {
        closeDiagramModal();
      }
    }

    document.addEventListener('keydown', function(e) {
      if (e.key === 'Escape') {
        closeDiagramModal();
        closeDbDetailModal();
      }
    });

    // Smart Paper Link Opener:
    // When viewed locally (file:// or localhost) and downloaded paper exists, opens local PDF/XML.
    // When served via GitHub Pages (where copyrighted papers are excluded via .gitignore),
    // automatically and gracefully redirects to the official publisher DOI URL.
    function smartOpenPaper(event, localUrl, doiUrl) {
      if (event) event.preventDefault();
      
      const isLocal = window.location.protocol === 'file:' || 
                      window.location.hostname === 'localhost' || 
                      window.location.hostname === '127.0.0.1' || 
                      window.location.hostname === '';
                      
      if (isLocal) {
        window.open(localUrl, '_blank');
      } else {
        window.open(doiUrl, '_blank');
      }
    }

    // Top View Mode Switcher (Summary vs Full Report vs DB Explorer vs BR-SAScore Theory)
    function switchViewTab(tabName, opts) {
      const btnIntro = document.getElementById('btn-tab-intro');
      const btnStoryline = document.getElementById('btn-tab-storyline');
      const btnTechnical = document.getElementById('btn-tab-technical');
      const btnDatabase = document.getElementById('btn-tab-database');
      const btnBrtheory = document.getElementById('btn-tab-brtheory');
      const btnGlossary = document.getElementById('btn-tab-glossary');
      const btnVerification = document.getElementById('btn-tab-verification');
      const contentIntro = document.getElementById('tab-content-intro');
      const contentStoryline = document.getElementById('tab-content-storyline');
      const contentTechnical = document.getElementById('tab-content-technical');
      const contentDatabase = document.getElementById('tab-content-database');
      const contentBrtheory = document.getElementById('tab-content-brtheory');
      const contentGlossary = document.getElementById('tab-content-glossary');
      const contentVerification = document.getElementById('tab-content-verification');
      const tocStoryline = document.getElementById('toc-storyline');
      const tocTechnical = document.getElementById('toc-technical');
      const tocDatabase = document.getElementById('toc-database');
      const tocBrtheory = document.getElementById('toc-brtheory');
      const tocSidebar = document.querySelector('.toc-sidebar');

      if (!btnStoryline || !btnTechnical) return;

      [btnIntro, btnStoryline, btnTechnical, btnDatabase, btnBrtheory, btnGlossary, btnVerification].forEach(b => b && b.classList.remove('active'));
      [contentIntro, contentStoryline, contentTechnical, contentDatabase, contentBrtheory, contentGlossary, contentVerification].forEach(c => { if (c) c.style.display = 'none'; });
      [tocStoryline, tocTechnical, tocBrtheory].forEach(t => { if (t) t.style.display = 'none'; });
      if (tocSidebar) tocSidebar.style.display = (tabName === 'glossary' || tabName === 'intro' || tabName === 'database' || tabName === 'verification') ? 'none' : '';

      const skipScrollTop = opts && opts.skipScrollTop;

      if (tabName === 'intro') {
        if (btnIntro) btnIntro.classList.add('active');
        if (contentIntro) contentIntro.style.display = 'block';
        if (!skipScrollTop) window.scrollTo(0, 0);
        requestAnimationFrame(() => {
          renderKatexInContainer('intro', contentIntro);
        });
      } else if (tabName === 'storyline') {
        btnStoryline.classList.add('active');
        if (contentStoryline) contentStoryline.style.display = 'block';
        if (tocStoryline) tocStoryline.style.display = 'block';
        if (!skipScrollTop) window.scrollTo(0, 0);
        requestAnimationFrame(() => {
          renderKatexInContainer('storyline', contentStoryline);
          renderMermaidInStorylineTab();
        });
      } else if (tabName === 'database') {
        if (btnDatabase) btnDatabase.classList.add('active');
        if (contentDatabase) contentDatabase.style.display = 'block';
        if (!skipScrollTop) window.scrollTo(0, 0);
        requestAnimationFrame(() => {
          ensureDbExplorerInitialized();
        });
      } else if (tabName === 'brtheory') {
        if (btnBrtheory) btnBrtheory.classList.add('active');
        if (contentBrtheory) contentBrtheory.style.display = 'block';
        if (tocBrtheory) tocBrtheory.style.display = 'block';
        if (!skipScrollTop) window.scrollTo(0, 0);
        requestAnimationFrame(() => {
          renderKatexInContainer('brtheory', contentBrtheory);
          renderMermaidInBrtheoryTab();
        });
      } else if (tabName === 'glossary') {
        if (btnGlossary) btnGlossary.classList.add('active');
        if (contentGlossary) contentGlossary.style.display = 'block';
        if (!skipScrollTop) window.scrollTo(0, 0);
        requestAnimationFrame(() => {
          renderKatexInContainer('glossary', contentGlossary);
        });
      } else if (tabName === 'verification') {
        if (btnVerification) btnVerification.classList.add('active');
        if (contentVerification) contentVerification.style.display = 'block';
        if (!skipScrollTop) window.scrollTo(0, 0);
        initVerifyTab();
      } else {
        btnTechnical.classList.add('active');
        if (contentTechnical) contentTechnical.style.display = 'block';
        if (tocTechnical) tocTechnical.style.display = 'block';
        if (!skipScrollTop) window.scrollTo(0, 0);
        requestAnimationFrame(() => {
          renderKatexInContainer('technical', contentTechnical);
          renderMermaidInTechnicalTab();
        });
      }
    }

    // -------------------------------------------------------------------------
    // Domain DB Explorer & 2D Structure Viewer State & Renderer
    // -------------------------------------------------------------------------
    let dbPayloadCache = null;
    let dbExplorerReady = false;
    let dbCurrentCategory = 'all';
    let dbCurrentViewMode = 'cards';
    const DB_PAGE_SIZE = 36;
    let dbFilteredRecords = [];
    let dbRenderedCount = 0;
    let dbSearchTimer = null;
    let dbLoadMoreObserver = null;

    function scheduleDbExplorerRender() {
      if (dbSearchTimer) clearTimeout(dbSearchTimer);
      dbSearchTimer = setTimeout(() => {
        dbSearchTimer = null;
        renderDbExplorer();
      }, 90);
    }

    function escapeHtmlStr(str) {
      return String(str == null ? '' : str)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;')
        .replace(/'/g, '&#39;');
    }

    function loadDbPayload() {
      if (dbPayloadCache) return dbPayloadCache;
      const el = document.getElementById('db-explorer-payload');
      if (!el) return { svgs: {}, records: [] };
      try {
        dbPayloadCache = JSON.parse(el.textContent || '{}');
      } catch (err) {
        console.error('DB payload parse error:', err);
        dbPayloadCache = { svgs: {}, records: [] };
      }
      return dbPayloadCache;
    }

    let dbCurrentMetal = 'all';
    const CORE_METAL_ORDER = ['Si', 'Ti', 'Zr', 'Hf', 'Al', 'Ta', 'Nb', 'W', 'Mo', 'Ru', 'Rh', 'Ir', 'Ga', 'In', 'As', 'Co', 'Ni', 'Cu', 'Pt', 'Zn', 'Sn', 'Ge', 'La', 'Y'];
    const NON_OTHER_ELEMENTS = new Set([...CORE_METAL_ORDER, 'Organic', '-', '']);
    const METAL_ALIAS_MAP = {
      si: 'Si', silicon: 'Si', 실리콘: 'Si', 규소: 'Si',
      ti: 'Ti', titanium: 'Ti', 티타늄: 'Ti', 타이타늄: 'Ti', 티탄: 'Ti',
      zr: 'Zr', zirconium: 'Zr', 지르코늄: 'Zr',
      hf: 'Hf', hafnium: 'Hf', 하프늄: 'Hf',
      al: 'Al', aluminum: 'Al', aluminium: 'Al', 알루미늄: 'Al', 알루미나: 'Al',
      ta: 'Ta', tantalum: 'Ta', 탄탈럼: 'Ta', 탄탈륨: 'Ta', 탄탈: 'Ta',
      w: 'W', tungsten: 'W', 텅스텐: 'W',
      mo: 'Mo', molybdenum: 'Mo', 몰리브덴: 'Mo', 몰리브데넘: 'Mo',
      ru: 'Ru', ruthenium: 'Ru', 루테늄: 'Ru',
      rh: 'Rh', rhodium: 'Rh', 로듐: 'Rh',
      ir: 'Ir', iridium: 'Ir', 이리듐: 'Ir',
      co: 'Co', cobalt: 'Co', 코발트: 'Co',
      sn: 'Sn', tin: 'Sn', 주석: 'Sn',
      zn: 'Zn', zinc: 'Zn', 아연: 'Zn',
      cu: 'Cu', copper: 'Cu', 구리: 'Cu',
      ni: 'Ni', nickel: 'Ni', 니켈: 'Ni',
      fe: 'Fe', iron: 'Fe', 철: 'Fe',
      mn: 'Mn', manganese: 'Mn', 망간: 'Mn', 망가니즈: 'Mn',
      nb: 'Nb', niobium: 'Nb', 니오븀: 'Nb',
      v: 'V', vanadium: 'V', 바나듐: 'V',
      in: 'In', indium: 'In', 인듐: 'In',
      ga: 'Ga', gallium: 'Ga', 갈륨: 'Ga',
      as: 'As', arsenic: 'As', 비소: 'As',
      ge: 'Ge', germanium: 'Ge', 게르마늄: 'Ge', 저마늄: 'Ge',
      la: 'La', lanthanum: 'La', 란타넘: 'La', 란탄: 'La',
      y: 'Y', yttrium: 'Y', 이트륨: 'Y',
      sr: 'Sr', strontium: 'Sr', 스트론튬: 'Sr',
      pt: 'Pt', platinum: 'Pt', 백금: 'Pt',
      pd: 'Pd', palladium: 'Pd', 팔라듐: 'Pd',
      li: 'Li', lithium: 'Li', 리튬: 'Li',
      mg: 'Mg', magnesium: 'Mg', 마그네슘: 'Mg',
      ca: 'Ca', calcium: 'Ca', 칼슘: 'Ca',
      ag: 'Ag', silver: 'Ag', 은: 'Ag',
      yb: 'Yb', er: 'Er', sc: 'Sc', cd: 'Cd', gd: 'Gd', b: 'B', p: 'P', se: 'Se', te: 'Te', sb: 'Sb', pb: 'Pb', na: 'Na', k: 'K',
      organic: 'Organic', 유기: 'Organic', 유기물: 'Organic', 비금속: 'Organic'
    };

    function resolveMetalQuery(rawQuery) {
      let q = (rawQuery || '').trim().toLowerCase();
      if (!q) return null;
      const prefixMatch = q.match(/^(?:metal|elem|element|금속|원소)\s*:\s*(.+)$/);
      if (prefixMatch) {
        q = prefixMatch[1].trim();
      }
      return METAL_ALIAS_MAP[q] || null;
    }

    const PERIODIC_TABLE_LAYOUT = [
      // [Z, symbol, row, col] — Periods 1-6 + Lanthanides (Row 8)
      [1,'H',1,1],[2,'He',1,18],
      [3,'Li',2,1],[4,'Be',2,2],[5,'B',2,13],[6,'C',2,14],[7,'N',2,15],[8,'O',2,16],[9,'F',2,17],[10,'Ne',2,18],
      [11,'Na',3,1],[12,'Mg',3,2],[13,'Al',3,13],[14,'Si',3,14],[15,'P',3,15],[16,'S',3,16],[17,'Cl',3,17],[18,'Ar',3,18],
      [19,'K',4,1],[20,'Ca',4,2],[21,'Sc',4,3],[22,'Ti',4,4],[23,'V',4,5],[24,'Cr',4,6],[25,'Mn',4,7],[26,'Fe',4,8],[27,'Co',4,9],[28,'Ni',4,10],[29,'Cu',4,11],[30,'Zn',4,12],[31,'Ga',4,13],[32,'Ge',4,14],[33,'As',4,15],[34,'Se',4,16],[35,'Br',4,17],[36,'Kr',4,18],
      [37,'Rb',5,1],[38,'Sr',5,2],[39,'Y',5,3],[40,'Zr',5,4],[41,'Nb',5,5],[42,'Mo',5,6],[43,'Tc',5,7],[44,'Ru',5,8],[45,'Rh',5,9],[46,'Pd',5,10],[47,'Ag',5,11],[48,'Cd',5,12],[49,'In',5,13],[50,'Sn',5,14],[51,'Sb',5,15],[52,'Te',5,16],[53,'I',5,17],[54,'Xe',5,18],
      [55,'Cs',6,1],[56,'Ba',6,2],[72,'Hf',6,4],[73,'Ta',6,5],[74,'W',6,6],[75,'Re',6,7],[76,'Os',6,8],[77,'Ir',6,9],[78,'Pt',6,10],[79,'Au',6,11],[80,'Hg',6,12],[81,'Tl',6,13],[82,'Pb',6,14],[83,'Bi',6,15],[84,'Po',6,16],[85,'At',6,17],[86,'Rn',6,18],
      // Lanthanides (row 8, cols 4-18)
      [57,'La',8,4],[58,'Ce',8,5],[59,'Pr',8,6],[60,'Nd',8,7],[61,'Pm',8,8],[62,'Sm',8,9],[63,'Eu',8,10],[64,'Gd',8,11],[65,'Tb',8,12],[66,'Dy',8,13],[67,'Ho',8,14],[68,'Er',8,15],[69,'Tm',8,16],[70,'Yb',8,17],[71,'Lu',8,18]
    ];

    let dbElemCountsCache = {};

    function getPeriodicHeatmapStyle(count, maxCount) {
      if (!count || count <= 0) return '';
      const t = maxCount > 1 ? Math.log(count) / Math.log(maxCount) : 1.0;
      const stops = [
        [0.00, [224, 242, 254]], // #e0f2fe
        [0.35, [125, 211, 252]], // #7dd3fc
        [0.68, [2, 132, 199]],   // #0284c7
        [1.00, [30, 58, 138]],   // #1e3a8a
      ];
      let c0 = stops[0], c1 = stops[stops.length - 1];
      for (let i = 0; i < stops.length - 1; i++) {
        if (t >= stops[i][0] && t <= stops[i + 1][0]) {
          c0 = stops[i];
          c1 = stops[i + 1];
          break;
        }
      }
      const f = (c1[0] === c0[0]) ? 0 : (t - c0[0]) / (c1[0] - c0[0]);
      const r = Math.round(c0[1][0] + (c1[1][0] - c0[1][0]) * f);
      const g = Math.round(c0[1][1] + (c1[1][1] - c0[1][1]) * f);
      const b = Math.round(c0[1][2] + (c1[1][2] - c0[1][2]) * f);
      const textCol = (t >= 0.48) ? '#ffffff' : '#0c4a6e';
      const brdCol = `rgb(${Math.max(0, r - 28)}, ${Math.max(0, g - 28)}, ${Math.max(0, b - 28)})`;
      return `background:rgb(${r},${g},${b}); color:${textCol}; border-color:${brdCol};`;
    }

    function renderPeriodicTableFilter(elemCounts, totalCount) {
      const grid = document.getElementById('dbPeriodicGrid');
      if (!grid) return;
      const ptElems = new Set(PERIODIC_TABLE_LAYOUT.map(e => e[1]));
      let maxCnt = 1;
      Object.entries(elemCounts).forEach(([k, v]) => {
        if (ptElems.has(k) && v > maxCnt) maxCnt = v;
      });
      const maxLbl = document.getElementById('dbPtMaxLabel');
      if (maxLbl) maxLbl.textContent = `${maxCnt}건`;
      const allBtn = document.getElementById('dbPtAllBtn');
      if (allBtn) allBtn.textContent = `전체 원소 (${totalCount})`;
      const orgBtn = document.getElementById('dbPtOrgBtn');
      if (orgBtn && elemCounts['Organic']) {
        orgBtn.style.display = '';
        orgBtn.textContent = `유기물 (${elemCounts['Organic']})`;
      }

      const cells = [];
      PERIODIC_TABLE_LAYOUT.forEach(([z, sym, row, col]) => {
        const cnt = elemCounts[sym] || 0;
        if (cnt > 0) {
          const hmStyle = getPeriodicHeatmapStyle(cnt, maxCnt);
          cells.push(
            `<button type="button" class="pt-cell has-data" data-dbmetal="${sym}" style="grid-row:${row}; grid-column:${col}; ${hmStyle}" title="${sym} (Z=${z}): ${cnt}건 — 클릭하여 필터링" onclick="toggleDbPeriodicElement('${sym}')">` +
              `<span class="pt-z">${z}</span><span class="pt-sym">${sym}</span><span class="pt-cnt">${cnt}</span>` +
            `</button>`
          );
        } else {
          cells.push(
            `<div class="pt-cell empty-data" style="grid-row:${row}; grid-column:${col};" aria-disabled="true">` +
              `<span class="pt-z">${z}</span><span class="pt-sym">${sym}</span>` +
            `</div>`
          );
        }
      });
      // Lanthanide placeholder at row 6, col 3 and label at row 8, cols 1-3
      cells.push(`<div class="pt-cell empty-data" style="grid-row:6; grid-column:3; font-size:0.62rem; font-weight:700;">57-71</div>`);
      cells.push(`<div class="pt-lanth-gap" style="grid-row:7;"></div>`);
      cells.push(`<div class="pt-lanth-lbl" style="grid-row:8;">Lanthanides &rarr;</div>`);
      grid.innerHTML = cells.join('');
    }

    function toggleDbPeriodicElement(sym) {
      setDbMetalFilter(dbCurrentMetal === sym ? 'all' : sym);
    }

    function ensureDbExplorerInitialized() {
      if (!dbExplorerReady) {
        const data = loadDbPayload();
        const elemCounts = {};
        if (data.records) {
          data.records.forEach(r => {
            const e = (r.element || '-').trim();
            if (e && e !== '-') {
              elemCounts[e] = (elemCounts[e] || 0) + 1;
            }
          });
        }
        dbElemCountsCache = elemCounts;
        renderPeriodicTableFilter(elemCounts, (data.records || []).length);
        const allBtn = document.getElementById('dbCatAllBtn');
        const precBtn = document.getElementById('dbCatPrecBtn');
        const expBtn = document.getElementById('dbCatExpandedBtn');
        const patBtn = document.getElementById('dbCatPatBtn');
        const inhBtn = document.getElementById('dbCatInhBtn');
        const reagBtn = document.getElementById('dbCatReagBtn');
        const rxnBtn = document.getElementById('dbCatRxnBtn');
        const negBtn = document.getElementById('dbCatNegBtn');
        const sideBadge = document.getElementById('dbSideBadge');
        const sideAll = document.getElementById('dbSideAllLink');
        const sidePrec = document.getElementById('dbSidePrecLink');
        const sideExp = document.getElementById('dbSideExpandedLink');
        const sidePat = document.getElementById('dbSidePatLink');
        const sideInh = document.getElementById('dbSideInhLink');
        const sideReag = document.getElementById('dbSideReagLink');
        const sideRxn = document.getElementById('dbSideRxnLink');
        const sideNeg = document.getElementById('dbSideNegLink');
        const sideSvg = document.getElementById('dbSideSvgCount');
        const svgCnt = data.svgs ? Object.keys(data.svgs).length : 0;
        const recs = data.records || [];
        const cntOf = (cat) => recs.filter(r => cat === 'expanded_master' ? r.is_expanded_master : (r.db_type === cat || (r.alt_db_types || []).includes(cat))).length;
        const cAll = recs.length;
        const cPrec = cntOf('precursor');
        const cExp = cntOf('expanded_master');
        const cPat = cntOf('patent_expansion');
        const cInh = cntOf('inhibitor');
        const cReag = cntOf('reagent');
        const cRxn = cntOf('reaction');
        const cNeg = cntOf('negative');
        if (allBtn) allBtn.textContent = `전체 통합 (${cAll})`;
        if (precBtn) precBtn.textContent = `반도체 공정 프리커서 Baseline (${cPrec})`;
        if (expBtn) expBtn.textContent = `확장 ALD 프리커서 DB (${cExp})`;
        if (patBtn) patBtn.textContent = `특허·문헌 확장 DB (${cPat})`;
        if (inhBtn) inhBtn.textContent = `AS-ALD 인히비터 (${cInh})`;
        if (reagBtn) reagBtn.textContent = `출발물질·시약 (${cReag})`;
        if (rxnBtn) rxnBtn.textContent = `공인 합성 반응 (${cRxn})`;
        if (negBtn) negBtn.textContent = `Negative Controls (${cNeg})`;
        if (sideBadge) sideBadge.textContent = `${cAll}건`;
        if (sideAll) sideAll.textContent = `전체 통합 DB (${cAll})`;
        if (sidePrec) sidePrec.textContent = `상용 프리커서 Baseline (${cPrec})`;
        if (sideExp) sideExp.textContent = `확장 ALD 프리커서 DB (${cExp})`;
        if (sidePat) sidePat.textContent = `특허·문헌 확장 DB (${cPat})`;
        if (sideInh) sideInh.textContent = `AS-ALD 인히비터 (${cInh})`;
        if (sideReag) sideReag.textContent = `출발물질·시약 (${cReag})`;
        if (sideRxn) sideRxn.textContent = `공인 합성 반응 (${cRxn})`;
        if (sideNeg) sideNeg.textContent = `Negative Controls (${cNeg})`;
        if (sideSvg) sideSvg.innerHTML = `고유 <code>SMILES</code> ${svgCnt}종 및 반응식 12종 실시간 시각화 지원 (정준 중복 0건).`;
        const headerStats = document.getElementById('dbHeaderStats');
        if (headerStats) {
          const n2d = recs.filter(r => r.has_2d).length;
          headerStats.innerHTML = `${cAll} Canonical Unique Entities (38 Patent/Lit Integrated) &bull; ${n2d} with 2D SVG (${svgCnt} Unique SMILES)`;
        }
        dbExplorerReady = true;
      }
      renderDbExplorer();
    }

    function setDbMetalFilter(metalVal) {
      dbCurrentMetal = metalVal || 'all';
      document.querySelectorAll('#dbPeriodicWrap [data-dbmetal]').forEach(btn => {
        btn.classList.toggle('active', btn.getAttribute('data-dbmetal') === dbCurrentMetal);
      });
      const badge = document.getElementById('dbPeriodicActiveBadge');
      if (badge) {
        if (dbCurrentMetal === 'all') {
          badge.textContent = '선택: 전체 원소';
        } else {
          const cnt = dbElemCountsCache[dbCurrentMetal] || 0;
          badge.textContent = `선택: ${dbCurrentMetal} (${cnt}건)`;
        }
      }
      renderDbExplorer();
    }

    function setDbCategory(cat) {
      dbCurrentCategory = cat;
      document.querySelectorAll('#dbCategoryPills .db-pill').forEach(btn => {
        btn.classList.toggle('active', btn.getAttribute('data-dbcat') === cat);
      });
      document.querySelectorAll('#toc-database .toc-link').forEach(link => {
        link.classList.toggle('active', link.getAttribute('data-dbcat') === cat);
      });
      renderDbExplorer();
    }

    function setDbCategoryFromSidebar(cat, el) {
      setDbCategory(cat);
    }

    function setDbViewMode(mode) {
      dbCurrentViewMode = mode;
      const bCards = document.getElementById('btnDbViewCards');
      const bTable = document.getElementById('btnDbViewTable');
      if (bCards) bCards.classList.toggle('active', mode === 'cards');
      if (bTable) bTable.classList.toggle('active', mode === 'table');
      renderDbExplorer();
    }

    function resetDbFilters() {
      const inp = document.getElementById('dbSearchInput');
      const modeSel = document.getElementById('dbSearchMode');
      const chk2d = document.getElementById('dbOnly2dCheck');
      const sortSel = document.getElementById('dbSortSelect');
      if (inp) inp.value = '';
      if (modeSel) modeSel.value = 'smart';
      if (chk2d) chk2d.checked = false;
      if (sortSel) sortSel.value = 'default';
      dbCurrentMetal = 'all';
      document.querySelectorAll('#dbPeriodicWrap [data-dbmetal]').forEach(btn => {
        btn.classList.toggle('active', btn.getAttribute('data-dbmetal') === 'all');
      });
      const badge = document.getElementById('dbPeriodicActiveBadge');
      if (badge) badge.textContent = '선택: 전체 원소';
      setDbCategory('all');
    }

    function getDbBadgeStyle(dbType) {
      if (dbType === 'precursor') return 'badge-prec';
      if (dbType === 'expanded') return 'badge-ext';
      if (dbType === 'patent_expansion') return 'badge-phys';
      if (dbType === 'inhibitor') return 'badge-2g';
      if (dbType === 'reagent') return 'badge-1g';
      if (dbType === 'reaction') return 'badge-3g';
      if (dbType === 'negative') return 'badge-4g';
      return 'badge-ext';
    }

    function buildCoreSearchBlob(r) {
      const parts = [
        r.id, r.name, r.chemical_name, r.formula, r.cas, r.element, r.category,
        r.smiles, r.target_films, r.reaction_type, r.atomiclimits_name,
        r.target_surface, growthSurfaceStr(r), r.jurisdictions, r.patent_references,
        (r.aliases || []).join(' '), (r.secondary_ids || []).join(' '),
        (r.vendors || []).join(' '), (r.starting_materials || []).join(' ')
      ];
      return parts.filter(Boolean).join(' ').toLowerCase();
    }

    function growthSurfaceStr(r) {
      return r.growth_surface || '';
    }

    function buildRecordSearchBlob(r) {
      const parts = [
        r.id, r.name, r.chemical_name, r.formula, r.cas, r.element, r.category,
        r.smiles, r.target_films, r.process_type, r.reaction_type, r.byproduct,
        r.solvent, r.conditions, r.target_surface, r.growth_surface,
        r.passivation_mechanism, r.surface_bond, r.handling, r.atomiclimits_name,
        r.jurisdictions, r.patent_references, r.mechanism_note, r.diag_note,
        (r.aliases || []).join(' '), (r.secondary_ids || []).join(' '), (r.role_badges || []).join(' '),
        (r.vendors || []).join(' '), (r.starting_materials || []).join(' '),
        (r.atomiclimits_materials || []).join(' '), (r.atomiclimits_co_reactants || []).join(' '),
        (r.reactants || []).join(' '), (r.products || []).join(' ')
      ];
      return parts.filter(Boolean).join(' ').toLowerCase();
    }

    function renderReactionSchemeHtml(r, svgs, compact) {
      const rx = r.reactants || [];
      const px = r.products || [];
      let html = '';
      rx.forEach((s, idx) => {
        if (idx > 0) html += '<span class="db-rxn-op">+</span>';
        const svg = svgs[s] || `<span style="font-family:'JetBrains Mono';font-size:0.7rem;">${escapeHtmlStr(s)}</span>`;
        html += `<div class="db-rxn-mol" title="${escapeHtmlStr(s)}">${svg}</div>`;
      });
      html += '<span class="db-rxn-op">&rarr;</span>';
      px.forEach((s, idx) => {
        if (idx > 0) html += '<span class="db-rxn-op">+</span>';
        const svg = svgs[s] || `<span style="font-family:'JetBrains Mono';font-size:0.7rem;">${escapeHtmlStr(s)}</span>`;
        html += `<div class="db-rxn-mol" title="${escapeHtmlStr(s)}">${svg}</div>`;
      });
      if (compact && r.byproduct) {
        html += `<span style="font-family:'JetBrains Mono';font-size:0.72rem;color:#64748b;margin-left:4px;">+ ${escapeHtmlStr(r.byproduct)}</span>`;
      }
      return html;
    }

    function renderDbExplorer() {
      const data = loadDbPayload();
      const svgs = data.svgs || {};
      const allRecords = data.records || [];
      const rawQuery = ((document.getElementById('dbSearchInput') || {}).value || '').trim();
      const query = rawQuery.toLowerCase();
      const searchMode = ((document.getElementById('dbSearchMode') || {}).value || 'smart');
      const only2d = Boolean((document.getElementById('dbOnly2dCheck') || {}).checked);
      const sortMode = ((document.getElementById('dbSortSelect') || {}).value || 'default');
      const resolvedMetal = (searchMode === 'smart' || searchMode === 'metal') ? resolveMetalQuery(rawQuery) : null;

      let filtered = allRecords.filter(r => {
        if (dbCurrentCategory === 'expanded_master') {
          if (!r.is_expanded_master) return false;
        } else if (dbCurrentCategory !== 'all' && r.db_type !== dbCurrentCategory && !(r.alt_db_types || []).includes(dbCurrentCategory)) {
          return false;
        }
        const elem = (r.element || '-').trim();
        if (dbCurrentMetal === 'OTHER_METALS') {
          if (NON_OTHER_ELEMENTS.has(elem)) return false;
        } else if (dbCurrentMetal !== 'all' && elem !== dbCurrentMetal) {
          return false;
        }
        if (only2d && !r.has_2d) return false;
        if (query) {
          if (searchMode === 'metal') {
            if (resolvedMetal) return elem === resolvedMetal;
            return elem.toLowerCase() === query;
          }
          if (searchMode === 'smart' && resolvedMetal) {
            return elem === resolvedMetal;
          }
          if (searchMode === 'smart') {
            return buildCoreSearchBlob(r).includes(query);
          }
          return buildRecordSearchBlob(r).includes(query);
        }
        return true;
      });

      const hintEl = document.getElementById('dbSearchHint');
      if (hintEl) {
        if (query && resolvedMetal && (searchMode === 'smart' || searchMode === 'metal')) {
          hintEl.style.display = 'inline-block';
          hintEl.textContent = `메인 금속 '${resolvedMetal}' 기준 ${filtered.length}건 필터링`;
        } else {
          hintEl.style.display = 'none';
          hintEl.textContent = '';
        }
      }

      if (sortMode === 'mw_asc') {
        filtered.sort((a, b) => (a.mw ?? 1e9) - (b.mw ?? 1e9));
      } else if (sortMode === 'mw_desc') {
        filtered.sort((a, b) => (b.mw ?? -1) - (a.mw ?? -1));
      } else if (sortMode === 'bde_desc') {
        filtered.sort((a, b) => (b.bde_min ?? -1e9) - (a.bde_min ?? -1e9));
      } else if (sortMode === 'sa_asc') {
        filtered.sort((a, b) => (a.sascore ?? 1e9) - (b.sascore ?? 1e9));
      }

      const cntBadge = document.getElementById('dbResultCountBadge');
      if (cntBadge) cntBadge.textContent = `검색 결과: ${filtered.length}건`;
      const sideBadge = document.getElementById('dbSideBadge');
      if (sideBadge) sideBadge.textContent = `${filtered.length}건`;

      const emptyEl = document.getElementById('dbEmptyState');
      const gridEl = document.getElementById('dbCardGrid');
      const tblWrap = document.getElementById('dbTableWrap');
      const tblBody = document.getElementById('dbTableBody');
      const loadMoreBar = document.getElementById('dbLoadMoreBar');

      window._dbFilteredMap = {};
      filtered.forEach((r, i) => { window._dbFilteredMap[i] = r; });
      dbFilteredRecords = filtered;
      dbRenderedCount = 0;

      if (filtered.length === 0) {
        if (emptyEl) emptyEl.style.display = 'block';
        if (gridEl) { gridEl.style.display = 'none'; gridEl.innerHTML = ''; }
        if (tblWrap) tblWrap.style.display = 'none';
        if (tblBody) tblBody.innerHTML = '';
        if (loadMoreBar) loadMoreBar.style.display = 'none';
        return;
      }
      if (emptyEl) emptyEl.style.display = 'none';

      if (dbCurrentViewMode === 'cards') {
        if (tblWrap) tblWrap.style.display = 'none';
        if (tblBody) tblBody.innerHTML = '';
        if (gridEl) {
          gridEl.style.display = 'grid';
          gridEl.innerHTML = '';
        }
      } else {
        if (gridEl) { gridEl.style.display = 'none'; gridEl.innerHTML = ''; }
        if (tblWrap) tblWrap.style.display = 'block';
        if (tblBody) tblBody.innerHTML = '';
      }

      appendDbRenderBatch(DB_PAGE_SIZE);
      ensureDbLoadMoreObserver();
    }

    function ensureDbLoadMoreObserver() {
      if (dbLoadMoreObserver || typeof IntersectionObserver === 'undefined') return;
      const sentinel = document.getElementById('dbLoadMoreSentinel');
      if (!sentinel) return;
      dbLoadMoreObserver = new IntersectionObserver((entries) => {
        for (const entry of entries) {
          if (entry.isIntersecting && dbRenderedCount < dbFilteredRecords.length) {
            appendDbRenderBatch(DB_PAGE_SIZE);
          }
        }
      }, { rootMargin: '400px 0px', threshold: 0.01 });
      dbLoadMoreObserver.observe(sentinel);
    }

    function appendDbRenderBatch(batchSize) {
      const data = loadDbPayload();
      const svgs = data.svgs || {};
      const total = dbFilteredRecords.length;
      if (dbRenderedCount >= total) {
        updateDbLoadMoreUi();
        return;
      }

      const startIdx = dbRenderedCount;
      const endIdx = Math.min(startIdx + (batchSize || DB_PAGE_SIZE), total);
      const slice = dbFilteredRecords.slice(startIdx, endIdx);
      dbRenderedCount = endIdx;

      if (dbCurrentViewMode === 'cards') {
        const gridEl = document.getElementById('dbCardGrid');
        if (gridEl) {
          const html = slice.map((r, offset) => {
            const idx = startIdx + offset;
            const badgeCls = getDbBadgeStyle(r.db_type);
            let visualHtml = '';
            if (r.db_type === 'reaction') {
              visualHtml = `<div class="db-rxn-stage">${renderReactionSchemeHtml(r, svgs, true)}</div>`;
            } else if (r.has_2d && svgs[r.smiles]) {
              visualHtml = `<div class="db-svg-stage">${svgs[r.smiles]}</div>`;
            } else {
              const ligTokens = (r.name || '').match(/\(([^()]+)\)/g) || [];
              const ligChips = ligTokens.map(t => `<span style="background:#e0f2fe;color:#0369a1;border:1px solid #bae6fd;padding:2px 7px;border-radius:999px;font-size:0.7rem;font-weight:700;font-family:'JetBrains Mono',monospace;">${escapeHtmlStr(t.slice(1, -1))}</span>`).join('');
              visualHtml = `<div class="db-svg-stage" style="background:linear-gradient(135deg,#f8fafc 0%,#f0f9ff 100%); color:#475569; font-size:0.76rem; text-align:center; flex-direction:column; gap:7px; padding:10px;">
                <div style="display:flex; align-items:center; justify-content:center; gap:6px; flex-wrap:wrap;">
                  <span style="display:inline-flex;align-items:center;justify-content:center;width:30px;height:30px;border-radius:50%;background:#0284c7;color:#fff;font-weight:800;font-size:0.8rem;font-family:'JetBrains Mono',monospace;">${escapeHtmlStr(r.element || 'M')}</span>
                  ${ligChips || `<span style="font-family:'JetBrains Mono',monospace;font-weight:700;color:#0f172a;">${escapeHtmlStr(r.name)}</span>`}
                </div>
                <span style="font-size:0.7rem; color:#64748b;">배위 약어 표기 (이성질체 미특정 &bull; AtomicLimits 문헌)</span>
              </div>`;
            }

            const sapisTxt = (r.sapis_score != null) ? Number(r.sapis_score).toFixed(2) : '-';
            let ensSubHtml = '';
            if (r.ens_n_models >= 2) {
              const gStyle = r.ens_grade === 'High'
                ? 'background:#dcfce7;color:#15803d;border:1px solid #86efac;'
                : (r.ens_grade === 'Medium'
                  ? 'background:#fef3c7;color:#b45309;border:1px solid #fde68a;'
                  : 'background:#fee2e2;color:#b91c1c;border:1px solid #fca5a5;');
              const gShort = r.ens_grade === 'High' ? 'High' : (r.ens_grade === 'Medium' ? 'Med' : 'DFT!');
              ensSubHtml = `<div style="display:flex;align-items:center;gap:4px;margin-top:2px;flex-wrap:wrap;"><span style="font-family:'JetBrains Mono',monospace;font-size:0.65rem;color:#64748b;">[${Number(r.ens_min).toFixed(2)}–${Number(r.ens_max).toFixed(2)}]</span><span style="font-size:0.62rem;font-weight:700;padding:0 5px;border-radius:999px;${gStyle}" title="4-uMLIP UQ 신뢰도 (σ_BDE=${r.ens_sigma_bde} kcal/mol)">${gShort}</span></div>`;
            }
            let verdictBadge = '<span class="db-mini-pill" style="background:#f1f5f9;color:#475569;border-color:#cbd5e1;">참고 DB</span>';
            if (r.sapis_pass === true) {
              verdictBadge = '<span class="db-mini-pill" style="background:#dcfce7;color:#15803d;border-color:#86efac;font-weight:700;">합성 가능</span>';
            } else if (r.sapis_pass === false) {
              verdictBadge = '<span class="db-mini-pill" style="background:#fee2e2;color:#b91c1c;border-color:#fca5a5;font-weight:700;">합성 불가</span>';
            }
            const srcTxt = escapeHtmlStr(r.card_source || r.db_label || '-');
            const riskTxt = escapeHtmlStr(r.card_risk || '-');
            const riskColor = (r.sapis_pass === false) ? '#b91c1c' : ((r.sapis_pass === true) ? '#15803d' : '#475569');

            const uniformSummaryHtml = `<div style="display:grid; grid-template-columns:1fr 1fr 1fr; gap:6px; background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px; padding:8px 10px; font-size:0.74rem; margin-top:4px;">
              <div>
                <div style="color:#64748b; font-size:0.67rem;">${r.ens_n_models >= 2 ? 'SAPIS (4-uMLIP)' : 'SAPIS 점수'}</div>
                <div style="font-family:'JetBrains Mono', monospace; font-weight:800; font-size:0.86rem; color:#0f172a; margin-top:1px;">${sapisTxt}</div>
                ${ensSubHtml}
              </div>
              <div>
                <div style="color:#64748b; font-size:0.67rem;">합성 판정</div>
                <div style="margin-top:2px;">${verdictBadge}</div>
              </div>
              <div>
                <div style="color:#64748b; font-size:0.67rem;">출처</div>
                <div style="font-weight:600; color:#334155; margin-top:1px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis;" title="${srcTxt}">${srcTxt}</div>
              </div>
              <div style="grid-column:1 / -1; border-top:1px dashed #e2e8f0; padding-top:5px;">
                <div style="color:#64748b; font-size:0.67rem;">${r.sapis_pass === false ? '불가 원인 / 위험성' : '위험성'}</div>
                <div style="font-weight:600; color:${riskColor}; margin-top:1px; line-height:1.4; word-break:keep-all;" title="${riskTxt}">${riskTxt}</div>
              </div>
            </div>`;

            let subLine = escapeHtmlStr(r.chemical_name);
            if (r.aliases && r.aliases.length) {
              subLine += ` &bull; 별칭: ${escapeHtmlStr(r.aliases.slice(0, 2).join(', '))}`;
            }
            if (r.db_type === 'precursor' || r.db_type === 'expanded' || r.db_type === 'patent_expansion') {
              if (r.target_surface) {
                subLine += ` &bull; 차단: ${escapeHtmlStr(r.target_surface)}`;
              } else {
                subLine += ` &bull; ${escapeHtmlStr(r.target_films || '-')}`;
              }
            } else if (r.db_type === 'inhibitor') {
              subLine += ` &bull; 차단: ${escapeHtmlStr(r.target_surface || '-')}`;
            } else if (r.db_type === 'reagent') {
              subLine = `${escapeHtmlStr(r.category)} &bull; ${escapeHtmlStr(r.purity || '-')}`;
            } else if (r.db_type === 'reaction') {
              subLine = `${escapeHtmlStr(r.solvent)} &bull; ${escapeHtmlStr(r.conditions)}`;
            }

            return `<div class="db-card ${r.db_type === 'reaction' ? 'db-card-reaction' : ''}" onclick="openDbDetailModal(${idx})">
              <div style="display:flex; justify-content:space-between; align-items:center; gap:6px;">
                <span class="badge ${badgeCls}" style="font-size:0.68rem;">${escapeHtmlStr(r.db_label)}</span>
                <div style="display:flex; gap:5px; align-items:center;">
                  ${r.element && r.element !== '-' ? `<span class="badge badge-1g" style="font-size:0.68rem;">${escapeHtmlStr(r.element)}</span>` : ''}
                  <span style="font-family:'JetBrains Mono'; font-size:0.71rem; color:#64748b;">${escapeHtmlStr(r.id)}${(r.secondary_ids && r.secondary_ids.length) ? ' (+' + r.secondary_ids.length + ')' : ''}</span>
                </div>
              </div>
              ${visualHtml}
              <div>
                <div style="font-weight:700; font-size:0.95rem; color:#0f172a; line-height:1.35;">${escapeHtmlStr(r.name)}</div>
                <div style="font-size:0.78rem; color:#475569; margin-top:2px; line-height:1.4;">${subLine}</div>
                <div style="font-family:'JetBrains Mono'; font-size:0.73rem; color:#0369a1; margin-top:3px;">${escapeHtmlStr(r.formula)}${r.cas && r.cas !== '-' ? ' &bull; CAS ' + escapeHtmlStr(r.cas) : ''}</div>
              </div>
              ${uniformSummaryHtml}
              ${r.smiles ? `<div style="font-family:'JetBrains Mono'; font-size:0.68rem; color:#64748b; background:#f8fafc; padding:4px 7px; border-radius:4px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis;" title="${escapeHtmlStr(r.smiles)}">${escapeHtmlStr(r.smiles)}</div>` : ''}
            </div>`;
          }).join('');
          gridEl.insertAdjacentHTML('beforeend', html);
        }
      } else {
        const tblBody = document.getElementById('dbTableBody');
        if (tblBody) {
          const html = slice.map((r, offset) => {
            const idx = startIdx + offset;
            const badgeCls = getDbBadgeStyle(r.db_type);
            let thumbHtml = '<span style="font-size:0.7rem; color:#94a3b8;">No SMILES</span>';
            if (r.db_type === 'reaction') {
              const p0 = (r.products || [])[0];
              thumbHtml = (p0 && svgs[p0]) ? svgs[p0] : '<span style="font-size:0.7rem; color:#0284c7;">Reaction</span>';
            } else if (r.has_2d && svgs[r.smiles]) {
              thumbHtml = svgs[r.smiles];
            }
            const metricTxt = [];
            if (r.mw != null) metricTxt.push(`MW ${r.mw}`);
            if (r.sascore != null) metricTxt.push(`SA ${r.sascore}`);
            if (r.scscore != null) metricTxt.push(`SC ${r.scscore}`);
            if (r.bde_min != null) metricTxt.push(`BDE ${r.bde_min}`);
            if (r.vbur != null) metricTxt.push(`%Vbur ${r.vbur}`);
            if (r.de_dimer != null) metricTxt.push(`dE_dim ${r.de_dimer}`);
            if (r.target_films && r.target_films !== '-') metricTxt.push(`Film: ${r.target_films}`);
            if (r.target_surface) metricTxt.push(`Passiv: ${r.target_surface}`);
            if (r.dg_kcal_mol != null) metricTxt.push(`dG ${r.dg_kcal_mol} kcal/mol`);

            return `<tr style="cursor:pointer;" onclick="openDbDetailModal(${idx})">
              <td class="db-thumb-cell" style="text-align:center;">${thumbHtml}</td>
              <td>
                <strong>${escapeHtmlStr(r.name)}</strong><br>
                <span style="font-family:'JetBrains Mono'; font-size:0.7rem; color:#64748b;">${escapeHtmlStr(r.id)}</span><br>
                <span class="badge ${badgeCls}" style="font-size:0.65rem; margin-top:2px;">${escapeHtmlStr(r.db_label)}</span>
              </td>
              <td>
                <div style="font-weight:600; font-size:0.76rem; color:#0f172a;">${escapeHtmlStr(r.chemical_name)}</div>
                <div style="font-family:'JetBrains Mono'; font-size:0.72rem; color:#0369a1;">${escapeHtmlStr(r.formula)}</div>
                <div style="font-family:'JetBrains Mono'; font-size:0.69rem; color:#64748b;">CAS: ${escapeHtmlStr(r.cas)}</div>
              </td>
              <td style="text-align:center;"><span class="badge badge-1g">${escapeHtmlStr(r.element || '-')}</span></td>
              <td style="font-size:0.74rem; line-height:1.45; color:#334155;">${escapeHtmlStr(metricTxt.join(' | ') || r.category || '-')}</td>
              <td style="font-family:'JetBrains Mono'; font-size:0.69rem; color:#0f172a;"><code>${escapeHtmlStr(r.smiles || '-')}</code></td>
            </tr>`;
          }).join('');
          tblBody.insertAdjacentHTML('beforeend', html);
        }
      }
      updateDbLoadMoreUi();
    }

    function updateDbLoadMoreUi() {
      const bar = document.getElementById('dbLoadMoreBar');
      const status = document.getElementById('dbLoadMoreStatus');
      const btnMore = document.getElementById('btnDbLoadMore');
      const total = dbFilteredRecords.length;
      if (!bar) return;
      if (total === 0 || dbRenderedCount >= total) {
        bar.style.display = 'none';
        return;
      }
      const remaining = total - dbRenderedCount;
      bar.style.display = 'flex';
      if (status) {
        status.textContent = `${dbRenderedCount} / ${total}건 표시 중 (남은 ${remaining}건은 스크롤 시 자동 로드)`;
      }
      if (btnMore) {
        btnMore.textContent = `더 보기 (+${Math.min(DB_PAGE_SIZE, remaining)}건)`;
      }
    }

    function copyDbText(btn, text) {
      if (!navigator.clipboard) return;
      navigator.clipboard.writeText(text).then(() => {
        const orig = btn.textContent;
        btn.textContent = '복사 완료';
        setTimeout(() => { btn.textContent = orig; }, 1400);
      });
    }

    function openDbDetailModal(filteredIdx) {
      const r = window._dbFilteredMap ? window._dbFilteredMap[filteredIdx] : null;
      if (!r) return;
      const data = loadDbPayload();
      const svgs = data.svgs || {};
      const modal = document.getElementById('dbDetailModal');
      const titleEl = document.getElementById('dbModalTitle');
      const subEl = document.getElementById('dbModalSub');
      const bodyEl = document.getElementById('dbModalBody');
      if (!modal || !bodyEl) return;

      titleEl.textContent = `${r.name} (${r.id})`;
      subEl.textContent = `${r.db_label} | ${r.chemical_name} | Formula: ${r.formula}${r.cas && r.cas !== '-' ? ' | CAS: ' + r.cas : ''}`;

      let svgBlock = '';
      if (r.db_type === 'reaction') {
        svgBlock = `<div class="db-detail-svg-wrap">${renderReactionSchemeHtml(r, svgs, false)}</div>`;
      } else if (r.has_2d && svgs[r.smiles]) {
        svgBlock = `<div class="db-detail-svg-wrap">${svgs[r.smiles]}</div>`;
      } else {
        svgBlock = `<div class="db-detail-svg-wrap" style="min-height:100px; color:#64748b;">2D SMILES 미수록 (AtomicLimits 문헌 기반 레코드)</div>`;
      }

      const rows = [];
      const addRow = (k, v) => {
        if (v !== undefined && v !== null && v !== '' && v !== '-') {
          rows.push(`<tr><th style="width:28%; text-align:left; background:#f8fafc;">${escapeHtmlStr(k)}</th><td>${v}</td></tr>`);
        }
      };

      addRow('Element / Category', `${escapeHtmlStr(r.element)} &bull; ${escapeHtmlStr(r.category)}`);
      if (r.secondary_ids && r.secondary_ids.length) addRow('Unified Cross-DB IDs', `<code>${escapeHtmlStr(r.id)}</code> &larr; Merged: <code>${escapeHtmlStr(r.secondary_ids.join(', '))}</code>${r.aliases && r.aliases.length ? ' | Aliases: ' + escapeHtmlStr(r.aliases.join(', ')) : ''}`);
      if (r.mw != null) addRow('Molecular Weight', `${r.mw} g/mol`);
      if (r.coord_number != null) addRow('Coordination / Valence e-', `CN<sub>eff</sub> = ${r.coord_number} &bull; N<sub>ve</sub> = ${r.valence_electrons}`);
      if (r.target_films) addRow('Target Films / Process', `${escapeHtmlStr(r.target_films)} (${escapeHtmlStr(r.process_type || '-')})`);
      if (r.jurisdictions || r.patent_references) addRow('Patent / Literature Expansion', `${r.jurisdictions ? '<strong>[' + escapeHtmlStr(r.jurisdictions) + ']</strong> ' : ''}${escapeHtmlStr(r.patent_references || '')}`);
      if (r.mechanism_note) addRow('Deposition / Chemical Note', escapeHtmlStr(r.mechanism_note));
      if (r.diag_note) addRow('SAPIS Boundary Diagnostic', escapeHtmlStr(r.diag_note));
      if (r.vapor_pressure) addRow('Vapor Pressure / T_decomp', `${escapeHtmlStr(r.vapor_pressure)} &bull; >${escapeHtmlStr(r.t_decomp_c)} &deg;C`);
      if (r.starting_materials && r.starting_materials.length) addRow('Starting Materials & Route', `${escapeHtmlStr(r.starting_materials.join(' + '))} (${escapeHtmlStr(r.reaction_type)}, ${escapeHtmlStr(r.yield_pct)}% yield)`);
      if (r.vendors && r.vendors.length) addRow('Commercial Vendors', escapeHtmlStr(r.vendors.join(', ')));
      if (r.target_surface) addRow('AS-ALD Selectivity', `Passivation: <strong>${escapeHtmlStr(r.target_surface)}</strong> | Growth: ${escapeHtmlStr(r.growth_surface)} (${escapeHtmlStr(r.selectivity_type)})`);
      if (r.passivation_mechanism) addRow('Mechanism & Surface Bond', `${escapeHtmlStr(r.passivation_mechanism)} (${escapeHtmlStr(r.surface_bond)})`);
      if (r.removal_method) addRow('Removal Method', escapeHtmlStr(r.removal_method));
      if (r.purity) addRow('Purity / Handling', `${escapeHtmlStr(r.purity)} &bull; ${escapeHtmlStr(r.handling)}`);
      if (r.solvent) addRow('Reaction Conditions', `Solvent: ${escapeHtmlStr(r.solvent)} | Conditions: ${escapeHtmlStr(r.conditions)} | Byproduct: ${escapeHtmlStr(r.byproduct)} | dG: ${r.dg_kcal_mol} kcal/mol`);
      if (r.atomiclimits_process_count) addRow('AtomicLimits ALD Data', `Processes: ${r.atomiclimits_process_count} | Materials: ${escapeHtmlStr((r.atomiclimits_materials || []).join(', '))} | Co-reactants: ${escapeHtmlStr((r.atomiclimits_co_reactants || []).join(', '))} | DOIs: ${r.atomiclimits_doi_count}`);
      if (r.sascore != null) addRow('Phase 1 SA Scores', `SAScore: <strong>${r.sascore}</strong> | SCScore: <strong>${r.scscore}</strong> | SYBA: <strong>${r.syba}</strong> | BR-SAScore: <strong>${r.br_sascore}</strong> | RAscore: <strong>${r.rascore}</strong> | Ligand SA: <strong>${r.ligand_sa}</strong>`);
      if (r.bde_min != null) addRow('SevenNet-0 Physics (LBFGS)', `BDE_min: <strong>${r.bde_min} kcal/mol</strong> | %V_bur: <strong>${r.vbur}%</strong>${r.vbur_d2_7net != null ? ' (D2 ' + r.vbur_d2_7net + '%, ETKDG ' + r.vbur_etkdg + '%)' : ''} | dE_dimer: <strong>${r.de_dimer} kcal/mol</strong>${r.dg_dimer_373 != null ? ' (dG373 ' + r.dg_dimer_373 + ')' : ''} | e_coh: <strong>${r.e_coh} eV/atom</strong>${r.de_form_bond != null ? ' | dEf/vM: <strong>' + r.de_form_bond + ' kcal/mol</strong>' : ''}`);
      if (r.tmc_sa_phys != null) addRow('Phase 2 Analytical Verdict', `SAPIS (7net): <strong>${r.tmc_sa_phys}</strong>${r.tmc_sa_no3d != null ? ' (ETKDG No-3D: ' + r.tmc_sa_no3d + ')' : ''} | Flag: <code>${escapeHtmlStr(r.domain_flag)}</code> | Pass: <strong>${r.phase2_pass}</strong>`);
      if (r.ens_n_models >= 2) {
        const pmStr = Object.entries(r.ens_per_model || {}).map(([k, v]) => `${escapeHtmlStr(k)}: <strong>${Number(v).toFixed(2)}</strong>`).join(' | ');
        addRow('4-uMLIP Ensemble (UQ)', `Median SAPIS: <strong>${Number(r.ens_sapis).toFixed(2)}</strong> [${Number(r.ens_min).toFixed(2)} ~ ${Number(r.ens_max).toFixed(2)}] | Confidence: <strong>${escapeHtmlStr(r.ens_grade)}</strong> (&sigma;<sub>BDE</sub> = ${r.ens_sigma_bde} kcal/mol, IQR = ${r.ens_iqr_bde} kcal/mol)<br><span style="font-size:0.76rem;color:#475569;">${pmStr}</span>`);
      }
      if (r.doi) addRow('Reference DOI', `<a href="https://doi.org/${escapeHtmlStr(r.doi)}" target="_blank">${escapeHtmlStr(r.doi)}</a>`);

      const smiBlock = r.smiles
        ? `<div style="display:flex; align-items:center; gap:8px; background:#f8fafc; border:1px solid #cbd5e1; border-radius:8px; padding:10px 12px;">
             <code style="flex:1; font-size:0.78rem; word-break:break-all; color:#0f172a;">${escapeHtmlStr(r.smiles)}</code>
             <button type="button" class="btn-modal-tool" id="btnCopyDbSmiles">복사</button>
           </div>`
        : '';

      bodyEl.innerHTML = `${svgBlock}${smiBlock}<table style="font-size:0.82rem;"><tbody>${rows.join('')}</tbody></table>`;
      const copyBtn = document.getElementById('btnCopyDbSmiles');
      if (copyBtn) {
        copyBtn.onclick = () => copyDbText(copyBtn, r.smiles);
      }
      modal.classList.add('active');
      document.body.style.overflow = 'hidden';
    }

    function closeDbDetailModal() {
      const modal = document.getElementById('dbDetailModal');
      if (modal) {
        modal.classList.remove('active');
        document.body.style.overflow = '';
      }
    }

    function closeDbDetailModalOnOverlay(event) {
      if (event.target && event.target.id === 'dbDetailModal') {
        closeDbDetailModal();
      }
    }

    let mermaidRenderSeq = 0;
    let mermaidObserver = null;

    async function renderSingleMermaidElement(el, prefix) {
      if (typeof mermaid === 'undefined' || !el) return;
      if (el.getAttribute('data-processed') === 'true' || el.dataset.rendering === 'true') return;
      const src = (el.dataset.mermaidSrc || el.textContent || '').trim();
      if (!src) return;
      el.dataset.rendering = 'true';
      el.dataset.mermaidSrc = src;
      try {
        mermaidRenderSeq += 1;
        const pfx = el.dataset.mermaidPrefix || prefix || 'mmd';
        const uniqueSvgId = `mermaid-${pfx}-${mermaidRenderSeq}-${(el.id || 'diag').replace(/[^a-zA-Z0-9_-]/g, '')}`;
        const { svg, bindFunctions } = await mermaid.render(uniqueSvgId, src);
        el.innerHTML = svg;
        el.setAttribute('data-processed', 'true');
        if (typeof bindFunctions === 'function') bindFunctions(el);
      } catch (e) {
        console.warn('Mermaid render notice:', el.id, e);
      } finally {
        delete el.dataset.rendering;
      }
    }

    function ensureMermaidObserver() {
      if (mermaidObserver || typeof IntersectionObserver === 'undefined') return mermaidObserver;
      mermaidObserver = new IntersectionObserver((entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            const el = entry.target;
            mermaidObserver.unobserve(el);
            renderSingleMermaidElement(el, el.dataset.mermaidPrefix || 'lazy');
          }
        });
      }, { rootMargin: '600px 0px', threshold: 0.01 });
      return mermaidObserver;
    }

    async function renderMermaidInContainer(selector, prefix) {
      if (typeof mermaid === 'undefined') return;
      const container = document.querySelector(selector);
      if (!container || container.offsetParent === null) return;
      const nodes = Array.from(container.querySelectorAll('.mermaid:not([data-processed="true"])'));
      if (nodes.length === 0) return;

      const obs = ensureMermaidObserver();
      if (!obs) {
        for (const el of nodes) {
          await renderSingleMermaidElement(el, prefix);
        }
        return;
      }
      nodes.forEach((el) => {
        el.dataset.mermaidPrefix = prefix;
        if (!el.dataset.mermaidObserved) {
          el.dataset.mermaidObserved = 'true';
          obs.observe(el);
        }
      });
    }

    function renderMermaidInStorylineTab() {
      renderMermaidInContainer('#tab-content-storyline', 'sum');
    }

    function renderMermaidInTechnicalTab() {
      renderMermaidInContainer('#tab-content-technical', 'tech');
    }

    function renderMermaidInBrtheoryTab() {
      renderMermaidInContainer('#tab-content-brtheory', 'brt');
    }

    function jumpToTechnicalSection(sectionId) {
      switchViewTab('technical', { skipScrollTop: true });
      requestAnimationFrame(() => {
        const el = document.getElementById(sectionId);
        if (el) el.scrollIntoView({ behavior: 'auto' });
      });
    }

    // ScrollSpy for Left Sticky TOC (Storyline, Technical & BR-SAScore Theory) - rAF throttled
    let glossaryActiveCat = 'all';
    function setGlossaryCategory(cat, btnEl) {
      glossaryActiveCat = cat || 'all';
      const wrap = document.getElementById('glossaryCatPills');
      if (wrap && btnEl) {
        wrap.querySelectorAll('button').forEach(b => b.classList.remove('active'));
        btnEl.classList.add('active');
      }
      filterGlossary();
    }
    function resetGlossaryFilter() {
      const inp = document.getElementById('glossarySearchInput');
      if (inp) inp.value = '';
      glossaryActiveCat = 'all';
      const wrap = document.getElementById('glossaryCatPills');
      if (wrap) {
        wrap.querySelectorAll('button').forEach((b, idx) => b.classList.toggle('active', idx === 0));
      }
      filterGlossary();
    }
    function filterGlossary() {
      const inp = document.getElementById('glossarySearchInput');
      const q = (inp ? inp.value : '').trim().toLowerCase();
      const cards = document.querySelectorAll('#glossaryCardGrid .glossary-card');
      let shown = 0;
      cards.forEach(c => {
        const cat = c.getAttribute('data-gcat') || '';
        const txt = c.getAttribute('data-gtext') || '';
        const okCat = (glossaryActiveCat === 'all') || (cat === glossaryActiveCat);
        const okTxt = !q || txt.indexOf(q) !== -1;
        const vis = okCat && okTxt;
        c.style.display = vis ? 'flex' : 'none';
        if (vis) shown += 1;
      });
      const badge = document.getElementById('glossaryCountBadge');
      if (badge) badge.textContent = '표시 ' + shown + '개';
      const empty = document.getElementById('glossaryEmptyState');
      if (empty) empty.style.display = shown === 0 ? 'block' : 'none';
    }
    function initScrollSpy() {
      const progressBadge = document.getElementById('tocProgress');
      const brtProgressBadge = document.getElementById('brtTocProgress');
      let scrollSpyTicking = false;
      let lastActiveId = null;
      let lastBadgeText = '';
      let lastBrtBadgeText = '';

      function updateScrollSpy() {
        scrollSpyTicking = false;
        const dbContainer = document.getElementById('tab-content-database');
        if (dbContainer && dbContainer.style.display !== 'none') return;

        const scrollPos = window.scrollY + 140;
        const storylineContainer = document.getElementById('tab-content-storyline');
        const brtheoryContainer = document.getElementById('tab-content-brtheory');
        const isStoryline = storylineContainer && storylineContainer.style.display !== 'none';
        const isBrtheory = brtheoryContainer && brtheoryContainer.style.display !== 'none';

        const activeContainer = isStoryline
          ? storylineContainer
          : (isBrtheory ? brtheoryContainer : document.getElementById('tab-content-technical'));
        if (!activeContainer) return;

        const sections = activeContainer.querySelectorAll('section[id]');
        const activeToc = isStoryline
          ? document.getElementById('toc-storyline')
          : (isBrtheory ? document.getElementById('toc-brtheory') : document.getElementById('toc-technical'));
        if (!activeToc || sections.length === 0) return;

        let currentIndex = 0;
        for (let idx = 0; idx < sections.length; idx++) {
          if (sections[idx].offsetTop <= scrollPos) {
            currentIndex = idx;
          } else {
            break;
          }
        }

        const activeId = sections[currentIndex] ? sections[currentIndex].getAttribute('id') : null;
        if (activeId && activeId !== lastActiveId) {
          lastActiveId = activeId;
          const tocLinks = activeToc.querySelectorAll('.toc-link');
          for (let i = 0; i < tocLinks.length; i++) {
            const link = tocLinks[i];
            const isMatch = link.getAttribute('data-target') === activeId;
            link.classList.toggle('active', isMatch);
          }
        }

        if (isBrtheory && brtProgressBadge) {
          const nextBrt = 'Sec ' + (currentIndex + 1) + '/' + sections.length;
          if (nextBrt !== lastBrtBadgeText) {
            lastBrtBadgeText = nextBrt;
            brtProgressBadge.textContent = nextBrt;
          }
        } else if (progressBadge && !isStoryline && !isBrtheory) {
          const nextBadge = 'Sec ' + (currentIndex + 1) + '/' + sections.length;
          if (nextBadge !== lastBadgeText) {
            lastBadgeText = nextBadge;
            progressBadge.textContent = nextBadge;
          }
        }
      }

      function onScroll() {
        if (!scrollSpyTicking) {
          scrollSpyTicking = true;
          requestAnimationFrame(updateScrollSpy);
        }
      }

      window.addEventListener('scroll', onScroll, { passive: true });
      updateScrollSpy();
    }

    (function wrapUncontainedTables() {
      document.addEventListener('DOMContentLoaded', function() {
        document.querySelectorAll('.main-content table').forEach(function(tbl) {
          var p = tbl.parentElement;
          if (!p) return;
          var st = p.getAttribute('style') || '';
          if (p.classList.contains('table-responsive') || p.classList.contains('mobile-table-scroll') || st.indexOf('overflow-x') !== -1) {
            return;
          }
          var wrap = document.createElement('div');
          wrap.className = 'mobile-table-scroll';
          p.insertBefore(wrap, tbl);
          wrap.appendChild(tbl);
        });
      });
    })();

    // =========================================================================
    // 304-Entity Dataset Provenance Verification Tab (localStorage + JSON I/O)
    // =========================================================================
    const VERIFY_STORAGE_KEY = 'sapis_source_verify_v1';
    let verifyTabInitialized = false;
    let verifyActiveCat = 'all';
    let verifyActiveStatus = 'all';
    let verifyState = { checked: {}, notes: {} };

    function loadVerifyState() {
      try {
        const raw = localStorage.getItem(VERIFY_STORAGE_KEY);
        if (raw) {
          const parsed = JSON.parse(raw);
          verifyState = {
            checked: (parsed && typeof parsed.checked === 'object') ? parsed.checked : {},
            notes: (parsed && typeof parsed.notes === 'object') ? parsed.notes : {},
          };
        }
      } catch (e) {
        verifyState = { checked: {}, notes: {} };
      }
    }

    function persistVerifyState() {
      try {
        localStorage.setItem(VERIFY_STORAGE_KEY, JSON.stringify(verifyState));
      } catch (e) {}
    }

    function updateVerifyProgressUi() {
      const rows = document.querySelectorAll('#verifyTableBody tr.ver-row');
      const total = rows.length;
      let done = 0;
      rows.forEach(tr => {
        const vid = tr.getAttribute('data-vid');
        const isDone = !!verifyState.checked[vid];
        if (isDone) done++;
        tr.style.background = isDone ? '#f0fdf4' : '';
      });
      const pct = total > 0 ? ((done * 100.0) / total).toFixed(1) : '0.0';
      const badge = document.getElementById('verifyProgressBadge');
      if (badge) badge.textContent = `확인 완료: ${done} / ${total} (${pct}%)`;
      const bar = document.getElementById('verifyProgressBar');
      if (bar) bar.style.width = `${pct}%`;
    }

    function initVerifyTab() {
      if (!verifyTabInitialized) {
        loadVerifyState();
        document.querySelectorAll('#verifyTableBody .ver-chk').forEach(chk => {
          const vid = chk.getAttribute('data-vid');
          chk.checked = !!verifyState.checked[vid];
        });
        document.querySelectorAll('#verifyTableBody .ver-note').forEach(inp => {
          const vid = inp.getAttribute('data-vid');
          inp.value = verifyState.notes[vid] || '';
        });
        verifyTabInitialized = true;
      }
      updateVerifyProgressUi();
      filterVerifyTable();
    }

    function toggleVerifyCheck(vid, isChecked) {
      if (isChecked) {
        verifyState.checked[vid] = true;
      } else {
        delete verifyState.checked[vid];
      }
      persistVerifyState();
      updateVerifyProgressUi();
      if (verifyActiveStatus !== 'all') filterVerifyTable();
    }

    function saveVerifyNote(vid, val) {
      const trimmed = (val || '').trim();
      if (trimmed) {
        verifyState.notes[vid] = trimmed;
      } else {
        delete verifyState.notes[vid];
      }
      persistVerifyState();
    }

    function setVerifyCategory(cat, btnEl) {
      verifyActiveCat = cat || 'all';
      const wrap = document.getElementById('verifyCatPills');
      if (wrap) wrap.querySelectorAll('.db-pill').forEach(b => b.classList.remove('active'));
      if (btnEl) btnEl.classList.add('active');
      filterVerifyTable();
    }

    function setVerifyStatus(st, btnEl) {
      verifyActiveStatus = st || 'all';
      const wrap = document.getElementById('verifyStatusPills');
      if (wrap) wrap.querySelectorAll('.db-pill').forEach(b => b.classList.remove('active'));
      if (btnEl) btnEl.classList.add('active');
      filterVerifyTable();
    }

    function filterVerifyTable() {
      const inp = document.getElementById('verifySearchInput');
      const q = (inp ? inp.value : '').trim().toLowerCase();
      const rows = document.querySelectorAll('#verifyTableBody tr.ver-row');
      let shown = 0;
      rows.forEach(tr => {
        const vid = tr.getAttribute('data-vid');
        const cat = tr.getAttribute('data-vcat');
        const blob = tr.getAttribute('data-vtext') || '';
        const note = (verifyState.notes[vid] || '').toLowerCase();
        const isDone = !!verifyState.checked[vid];
        const okCat = (verifyActiveCat === 'all') || (cat === verifyActiveCat);
        const okStat = (verifyActiveStatus === 'all') || (verifyActiveStatus === 'checked' && isDone) || (verifyActiveStatus === 'unchecked' && !isDone);
        const okQ = !q || blob.includes(q) || note.includes(q);
        const vis = okCat && okStat && okQ;
        tr.style.display = vis ? '' : 'none';
        if (vis) shown++;
      });
      const shownBadge = document.getElementById('verifyShownBadge');
      if (shownBadge) shownBadge.textContent = `표시: ${shown}건`;
    }

    function exportVerifyState() {
      const payload = {
        exported_at: new Date().toISOString(),
        checked_count: Object.keys(verifyState.checked).length,
        checked: verifyState.checked,
        notes: verifyState.notes,
      };
      const blob = new Blob([JSON.stringify(payload, null, 2)], { type: 'application/json' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'sapis_304_verification_log.json';
      a.click();
      URL.revokeObjectURL(url);
    }

    function importVerifyState(evt) {
      const file = evt && evt.target && evt.target.files && evt.target.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = function(e) {
        try {
          const parsed = JSON.parse(e.target.result);
          verifyState = {
            checked: (parsed && typeof parsed.checked === 'object') ? parsed.checked : {},
            notes: (parsed && typeof parsed.notes === 'object') ? parsed.notes : {},
          };
          persistVerifyState();
          verifyTabInitialized = false;
          initVerifyTab();
        } catch (err) {
          alert('유효한 JSON 파일이 아닙니다.');
        }
      };
      reader.readAsText(file);
    }

    function resetVerifyChecks() {
      if (!confirm('모든 검증 체크 및 메모 기록을 초기화하시겠습니까?')) return;
      verifyState = { checked: {}, notes: {} };
      persistVerifyState();
      verifyTabInitialized = false;
      initVerifyTab();
    }
