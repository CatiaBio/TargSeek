#!/usr/bin/env python3
"""
Shared HTML shell for TargSeek reports.

Produces a single self-contained page (no external assets, so it opens from the
filesystem and survives being emailed) with a sticky sidebar for jumping between
sections, light/dark support, and print styles.

Colours follow a validated palette:
  - sequential blue for score magnitude
  - the fixed status palette for recommendations, always paired with a label and
    an icon so meaning never rests on colour alone
  - hairline rules and recessive greys rather than heavy borders

Used by integrate_comprehensive_epitope_analysis.py; kept separate so other
reports can reuse the same shell.
"""

import html
from datetime import datetime

# --- status palette (fixed, never themed) -----------------------------------
STATUS = {
    "good":     "#0aac23",
    "warning":  "#fbb315",
    "serious":  "#f18a59",
    "critical": "#dc3c3d",
}

# Recommendation tier -> (status role, filled?, icon name, short label)
RECOMMENDATION_TIERS = [
    ("Highly Recommended", "good", True, "check-double", "Highly recommended"),
    ("Recommended", "good", False, "check", "Recommended"),
    ("Moderate", "warning", False, "dash", "Moderate"),
    ("Low Priority", "critical", False, "cross", "Low priority"),
]

ICONS = {
    # 12x12 viewBox, currentColor, no fill-rule surprises
    "check-double": '<path d="M1 6.5 3.5 9 8 3.5" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/><path d="M6.5 8.6 8 10l3.2-4" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>',
    "check": '<path d="M2 6.5 5 9.5 10 3" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>',
    "dash": '<path d="M2.5 6h7" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"/>',
    "cross": '<path d="M3 3l6 6M9 3l-6 6" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"/>',
    "info": '<circle cx="6" cy="6" r="5" fill="none" stroke="currentColor" stroke-width="1.5"/><path d="M6 5.2v3.3M6 3.6v.1" stroke="currentColor" stroke-width="1.6" stroke-linecap="round"/>',
    "alert": '<path d="M6 1.5 11 10.5H1z" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linejoin="round"/><path d="M6 5v2.4M6 8.8v.1" stroke="currentColor" stroke-width="1.5" stroke-linecap="round"/>',
}


def esc(value):
    """HTML-escape any value."""
    return html.escape("" if value is None else str(value))


def icon(name, size=12):
    body = ICONS.get(name, ICONS["info"])
    return (f'<svg class="ic" width="{size}" height="{size}" viewBox="0 0 12 12" '
            f'aria-hidden="true" focusable="false">{body}</svg>')


# --- small components -------------------------------------------------------

def stat_tile(value, label, hint=None):
    """One figure with its label. Not a chart: a single number reads faster."""
    hint_html = f'<div class="tile-hint">{esc(hint)}</div>' if hint else ""
    return (f'<div class="tile"><div class="tile-value">{esc(value)}</div>'
            f'<div class="tile-label">{esc(label)}</div>{hint_html}</div>')


def stat_row(tiles):
    return f'<div class="tiles">{"".join(tiles)}</div>'


def score_num(value, decimals=2):
    """Plain number for a component score.

    Component columns are read by comparison down the column, which tabular
    figures already do well; a bar in every cell of every column makes the table
    wider than the screen and adds no information.
    """
    try:
        return f'<span class="scorenum">{float(value):.{decimals}f}</span>'
    except (TypeError, ValueError):
        return '<span class="muted">-</span>'


def score_meter(value, decimals=2):
    """Inline meter for a 0-1 score.

    The number is always printed next to the bar, so the value never depends on
    reading a length or a colour.
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return '<span class="muted">-</span>'
    pct = max(0.0, min(1.0, v)) * 100
    return (f'<span class="meter" role="img" aria-label="{v:.{decimals}f} of 1">'
            f'<span class="meter-track"><span class="meter-fill" style="width:{pct:.1f}%"></span></span>'
            f'<span class="meter-num">{v:.{decimals}f}</span></span>')


def recommendation_badge(text):
    """Status badge: colour + icon + words, never colour alone."""
    raw = "" if text is None else str(text)
    for needle, role, filled, icon_name, label in RECOMMENDATION_TIERS:
        if needle.lower() in raw.lower():
            cls = f"badge badge-{role}" + (" badge-filled" if filled else "")
            return f'<span class="{cls}">{icon(icon_name)}{esc(label)}</span>'
    return f'<span class="badge">{icon("info")}{esc(raw or "n/a")}</span>'


def note(kind, title, body_html):
    """Callout. kind: info | warning | critical | good."""
    icon_name = {"warning": "alert", "critical": "alert", "good": "check"}.get(kind, "info")
    return (f'<div class="note note-{kind}"><div class="note-head">{icon(icon_name, 14)}'
            f'<strong>{esc(title)}</strong></div><div class="note-body">{body_html}</div></div>')


def table(headers, rows, sortable=True, table_id=None, numeric_from=None):
    """Data table.

    headers: list of column names
    rows:    list of lists of pre-rendered HTML cells
    numeric_from: index at which columns become numeric (right-aligned)
    """
    tid = f' id="{table_id}"' if table_id else ""
    cls = "data-table" + (" sortable" if sortable else "")
    head_cells = []
    for i, h in enumerate(headers):
        align = ' class="num"' if numeric_from is not None and i >= numeric_from else ""
        if sortable:
            head_cells.append(
                f'<th{align}><button type="button" class="sort-btn" data-col="{i}">'
                f'{esc(h)}<span class="sort-ind" aria-hidden="true"></span></button></th>')
        else:
            head_cells.append(f'<th{align}>{esc(h)}</th>')
    body = []
    for row in rows:
        cells = []
        for i, cell in enumerate(row):
            align = ' class="num"' if numeric_from is not None and i >= numeric_from else ""
            cells.append(f"<td{align}>{cell}</td>")
        body.append(f"<tr>{''.join(cells)}</tr>")
    return (f'<div class="table-wrap"><table{tid} class="{cls}">'
            f'<thead><tr>{"".join(head_cells)}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table></div>')


def section(section_id, title, body_html, subtitle=None, level=2):
    sub = f'<p class="section-sub">{esc(subtitle)}</p>' if subtitle else ""
    return (f'<section id="{section_id}" class="section">'
            f'<h{level}>{esc(title)}</h{level}>{sub}{body_html}</section>')


# --- the page --------------------------------------------------------------

CSS = """
:root {
  color-scheme: light;
  --surface:        #ffffff;
  --plane:          #f9f9f9;
  --sidebar:        #ffffff;
  --text:           #191919;
  --text-2:         #666666;
  --hairline:       #e7e7e7;
  --border:         rgba(11,11,11,0.10);
  --accent:         #257de0;
  --accent-soft:    #cce3fd;
  --track:          #eef1f4;
  --shadow:         0 1px 2px rgba(11,11,11,0.04), 0 8px 24px rgba(11,11,11,0.04);
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --surface:      #191919;
    --plane:        #0b0b0b;
    --sidebar:      #141414;
    --text:         #efefef;
    --text-2:       #9c9c9c;
    --hairline:     #303030;
    --border:       rgba(255,255,255,0.10);
    --accent:       #338ced;
    --accent-soft:  #1b64be;
    --track:        #242424;
    --shadow:       none;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface:      #191919;
  --plane:        #0b0b0b;
  --sidebar:      #141414;
  --text:         #efefef;
  --text-2:       #9c9c9c;
  --hairline:     #303030;
  --border:       rgba(255,255,255,0.10);
  --accent:       #338ced;
  --accent-soft:  #1b64be;
  --track:        #242424;
  --shadow:       none;
}

* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body {
  margin: 0;
  background: var(--plane);
  color: var(--text);
  font: 15px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif;
  -webkit-font-smoothing: antialiased;
}
a { color: var(--accent); }

/* layout */
.layout { display: grid; grid-template-columns: 264px minmax(0,1fr); min-height: 100vh; }
.sidebar {
  position: sticky; top: 0; align-self: start; height: 100vh;
  background: var(--sidebar); border-right: 1px solid var(--hairline);
  display: flex; flex-direction: column; overflow: hidden;
}
.sb-head { padding: 22px 20px 16px; border-bottom: 1px solid var(--hairline); }
.sb-brand { font-size: 12px; letter-spacing: .09em; text-transform: uppercase; color: var(--text-2); font-weight: 600; }
.sb-title { font-size: 17px; font-weight: 650; margin: 7px 0 10px; line-height: 1.3; }
.sb-meta { font-size: 12.5px; color: var(--text-2); }
.sb-meta div { margin-top: 3px; }
.sb-nav { padding: 12px 10px 20px; overflow-y: auto; flex: 1; }
.sb-group { margin-top: 14px; }
.sb-group-label {
  font-size: 11px; letter-spacing: .08em; text-transform: uppercase;
  color: var(--text-2); font-weight: 600; padding: 0 10px 6px;
}
.sb-nav a {
  display: flex; align-items: center; gap: 8px; justify-content: space-between;
  padding: 7px 10px; border-radius: 7px; color: var(--text);
  text-decoration: none; font-size: 13.5px; border-left: 2px solid transparent;
}
.sb-nav a:hover { background: var(--plane); }
.sb-nav a.active { background: var(--plane); border-left-color: var(--accent); font-weight: 600; }
.sb-nav a .count { font-size: 11.5px; color: var(--text-2); font-variant-numeric: tabular-nums; }
.sb-foot { padding: 12px 16px; border-top: 1px solid var(--hairline); font-size: 12px; color: var(--text-2); }

.main { padding: 0 0 72px; min-width: 0; }
.topbar {
  position: sticky; top: 0; z-index: 5;
  display: flex; align-items: center; gap: 12px; flex-wrap: wrap;
  padding: 12px 32px; background: color-mix(in srgb, var(--plane) 88%, transparent);
  backdrop-filter: blur(8px); border-bottom: 1px solid var(--hairline);
}
.topbar h1 { font-size: 15px; font-weight: 650; margin: 0; flex: 1; min-width: 180px; }
.content { padding: 8px 32px 0; max-width: 1180px; }

/* sections */
.section {
  background: var(--surface); border: 1px solid var(--border); border-radius: 12px;
  padding: 22px 24px; margin: 18px 0; box-shadow: var(--shadow);
  scroll-margin-top: 72px;
}
.section > h2 { font-size: 17px; font-weight: 650; margin: 0 0 4px; letter-spacing: -0.01em; }
.section > h3 { font-size: 14.5px; font-weight: 650; margin: 20px 0 8px; }
.section-sub { color: var(--text-2); font-size: 13.5px; margin: 0 0 16px; }
.section p { margin: 10px 0; }

/* stat tiles */
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(146px, 1fr)); gap: 12px; margin: 4px 0 6px; }
.tile { padding: 14px 16px; border: 1px solid var(--hairline); border-radius: 10px; background: var(--plane); }
.tile-value { font-size: 30px; font-weight: 640; line-height: 1.05; letter-spacing: -0.02em; }
.tile-label { font-size: 12.5px; color: var(--text-2); margin-top: 5px; }
.tile-hint { font-size: 11.5px; color: var(--text-2); margin-top: 3px; }

/* tables */
.table-wrap { overflow-x: auto; margin: 12px 0 4px; border: 1px solid var(--hairline); border-radius: 10px; }
table.data-table { width: 100%; border-collapse: collapse; font-size: 13px; }
.data-table th, .data-table td { padding: 9px 10px; text-align: left; border-bottom: 1px solid var(--hairline); white-space: nowrap; }
.data-table thead th { position: sticky; top: 0; background: var(--surface); z-index: 1; font-weight: 600; font-size: 12px; color: var(--text-2); }
.data-table tbody tr:last-child td { border-bottom: 0; }
.data-table tbody tr:hover { background: var(--plane); }
.data-table td.num, .data-table th.num { text-align: right; font-variant-numeric: tabular-nums; }
.sort-btn { all: unset; cursor: pointer; display: inline-flex; align-items: center; gap: 5px; font: inherit; color: inherit; }
.sort-btn:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; border-radius: 3px; }
.sort-ind::after { content: "\\2195"; opacity: .35; font-size: 11px; }
th[aria-sort="ascending"] .sort-ind::after { content: "\\2191"; opacity: 1; }
th[aria-sort="descending"] .sort-ind::after { content: "\\2193"; opacity: 1; }
.scorenum { font-variant-numeric: tabular-nums; }
.seq { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12.5px; letter-spacing: .02em; }

/* meter: thin bar, square at the baseline, 4px rounded data-end */
.meter { display: inline-flex; align-items: center; gap: 8px; }
.meter-track { width: 46px; height: 6px; background: var(--track); border-radius: 3px; overflow: hidden; flex: none; }
.meter-fill { display: block; height: 100%; background: var(--accent); border-radius: 0 4px 4px 0; }
.meter-num { font-variant-numeric: tabular-nums; min-width: 2.6em; text-align: right; }

/* badges */
.badge { display: inline-flex; align-items: center; gap: 5px; padding: 3px 9px; border-radius: 999px; font-size: 12px; font-weight: 600; white-space: nowrap; border: 1px solid var(--border); color: var(--text); }
.badge .ic { flex: none; }
.badge-good { color: #006c10; border-color: color-mix(in srgb, #0aac23 40%, transparent); background: color-mix(in srgb, #0aac23 10%, transparent); }
.badge-warning { color: #7a5200; border-color: color-mix(in srgb, #fbb315 45%, transparent); background: color-mix(in srgb, #fbb315 14%, transparent); }
.badge-serious { color: #8a3b14; border-color: color-mix(in srgb, #f18a59 45%, transparent); background: color-mix(in srgb, #f18a59 14%, transparent); }
.badge-critical { color: #a51f20; border-color: color-mix(in srgb, #dc3c3d 40%, transparent); background: color-mix(in srgb, #dc3c3d 10%, transparent); }
.badge-filled { background: #0aac23; border-color: #0aac23; color: #ffffff; }
:root[data-theme="dark"] .badge-good, :root[data-theme="dark"] .badge-warning,
:root[data-theme="dark"] .badge-serious, :root[data-theme="dark"] .badge-critical { color: var(--text); }
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) .badge-good,
  :root:where(:not([data-theme="light"])) .badge-warning,
  :root:where(:not([data-theme="light"])) .badge-serious,
  :root:where(:not([data-theme="light"])) .badge-critical { color: var(--text); }
}

/* notes */
.note { border: 1px solid var(--border); border-left-width: 3px; border-radius: 8px; padding: 12px 14px; margin: 12px 0; background: var(--plane); }
.note-head { display: flex; align-items: center; gap: 7px; margin-bottom: 4px; }
.note-body { font-size: 13.5px; color: var(--text); }
.note-body ul { margin: 6px 0 0; padding-left: 20px; }
.note-info { border-left-color: var(--accent); }
.note-good { border-left-color: #0aac23; }
.note-warning { border-left-color: #fbb315; }
.note-critical { border-left-color: #dc3c3d; }

/* controls */
.ctrl { display: inline-flex; align-items: center; gap: 7px; }
.search {
  font: inherit; font-size: 13px; padding: 6px 11px; min-width: 210px;
  color: var(--text); background: var(--surface);
  border: 1px solid var(--hairline); border-radius: 8px;
}
.search:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
.btn {
  font: inherit; font-size: 13px; padding: 6px 11px; cursor: pointer;
  color: var(--text); background: var(--surface);
  border: 1px solid var(--hairline); border-radius: 8px;
}
.btn:hover { background: var(--plane); }
.btn:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }
.muted { color: var(--text-2); }
.dl { display: grid; grid-template-columns: max-content 1fr; gap: 4px 16px; font-size: 13.5px; margin: 8px 0 0; }
.dl dt { color: var(--text-2); }
.dl dd { margin: 0; }
.protein-head { display: flex; align-items: baseline; gap: 12px; flex-wrap: wrap; margin-top: 2px; }
.protein-head h3 { margin: 0; }
.chiprow { display: flex; gap: 8px; flex-wrap: wrap; font-size: 12.5px; color: var(--text-2); }
.chip { border: 1px solid var(--hairline); border-radius: 999px; padding: 2px 9px; }

/* narrow screens: sidebar becomes a strip on top */
@media (max-width: 900px) {
  .layout { grid-template-columns: 1fr; }
  .sidebar { position: static; height: auto; border-right: 0; border-bottom: 1px solid var(--hairline); }
  .sb-nav { max-height: 210px; }
  .content, .topbar { padding-left: 18px; padding-right: 18px; }
}

/* print */
@media print {
  .sidebar, .topbar { display: none !important; }
  .layout { display: block; }
  body { background: #fff; }
  .section { break-inside: avoid; box-shadow: none; border-color: #ddd; margin: 10px 0; }
  .data-table th, .data-table td { white-space: normal; }
  .table-wrap { overflow: visible; }
}
"""

JS = """
// scroll-spy: highlight the section currently in view
const links = Array.from(document.querySelectorAll('.sb-nav a[href^="#"]'));
const byId = new Map(links.map(a => [a.getAttribute('href').slice(1), a]));
const seen = new Set();
const spy = new IntersectionObserver((entries) => {
  entries.forEach(e => e.isIntersecting ? seen.add(e.target.id) : seen.delete(e.target.id));
  const order = Array.from(document.querySelectorAll('.section')).map(s => s.id);
  const current = order.find(id => seen.has(id));
  links.forEach(a => a.classList.remove('active'));
  if (current && byId.has(current)) byId.get(current).classList.add('active');
}, { rootMargin: '-72px 0px -65% 0px', threshold: 0 });
document.querySelectorAll('.section').forEach(s => spy.observe(s));

// theme toggle, remembered per browser
const themeBtn = document.getElementById('theme-toggle');
if (themeBtn) {
  const saved = localStorage.getItem('targseek-theme');
  if (saved) document.documentElement.dataset.theme = saved;
  themeBtn.addEventListener('click', () => {
    const dark = document.documentElement.dataset.theme === 'dark'
      || (!document.documentElement.dataset.theme
          && window.matchMedia('(prefers-color-scheme: dark)').matches);
    const next = dark ? 'light' : 'dark';
    document.documentElement.dataset.theme = next;
    localStorage.setItem('targseek-theme', next);
  });
}

// filter every table by text
const search = document.getElementById('table-search');
if (search) {
  search.addEventListener('input', () => {
    const q = search.value.trim().toLowerCase();
    document.querySelectorAll('table.data-table tbody').forEach(tb => {
      let shown = 0;
      tb.querySelectorAll('tr').forEach(tr => {
        const hit = !q || tr.textContent.toLowerCase().includes(q);
        tr.hidden = !hit;
        if (hit) shown++;
      });
      const empty = tb.parentElement.querySelector('.no-match');
      if (empty) empty.hidden = shown > 0;
    });
  });
}

// click a header to sort
document.querySelectorAll('table.sortable').forEach(tbl => {
  tbl.querySelectorAll('th .sort-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      const col = Number(btn.dataset.col);
      const th = btn.closest('th');
      const desc = th.getAttribute('aria-sort') !== 'descending';
      tbl.querySelectorAll('th').forEach(h => h.removeAttribute('aria-sort'));
      th.setAttribute('aria-sort', desc ? 'descending' : 'ascending');
      const body = tbl.tBodies[0];
      const rows = Array.from(body.rows);
      const val = (tr) => {
        const cell = tr.cells[col];
        if (!cell) return '';
        const n = parseFloat(cell.textContent.replace(/[^0-9.\\-]/g, ''));
        return isNaN(n) ? cell.textContent.trim().toLowerCase() : n;
      };
      rows.sort((a, b) => {
        const x = val(a), y = val(b);
        if (typeof x === 'number' && typeof y === 'number') return desc ? y - x : x - y;
        return desc ? String(y).localeCompare(String(x)) : String(x).localeCompare(String(y));
      });
      rows.forEach(r => body.appendChild(r));
    });
  });
});
"""


def render_page(title, subtitle, meta_items, nav_groups, sections_html,
                brand="TargSeek", footer_note=None, show_controls=True):
    """Assemble the full page.

    meta_items:  list of (label, value) shown in the sidebar head
    nav_groups:  list of (group_label, [(anchor_id, label, count_or_None), ...])
    sections_html: the already-rendered <section> blocks
    """
    meta_html = "".join(
        f"<div><span class='muted'>{esc(k)}:</span> {esc(v)}</div>" for k, v in meta_items)

    nav_html = []
    for group_label, items in nav_groups:
        links = []
        for anchor, label, count in items:
            badge = f"<span class='count'>{esc(count)}</span>" if count is not None else ""
            links.append(f"<a href='#{esc(anchor)}'><span>{esc(label)}</span>{badge}</a>")
        label_html = (f"<div class='sb-group-label'>{esc(group_label)}</div>"
                      if group_label else "")
        nav_html.append(f"<div class='sb-group'>{label_html}{''.join(links)}</div>")

    controls = ""
    if show_controls:
        controls = (
            '<div class="ctrl">'
            '<label class="sr-only" for="table-search">Filter tables</label>'
            '<input id="table-search" class="search" type="search" '
            'placeholder="Filter tables (gene, sequence, ...)">'
            '</div>'
            '<button id="theme-toggle" class="btn" type="button">Light / dark</button>'
            '<button class="btn" type="button" onclick="window.print()">Print</button>'
        )

    foot = footer_note or f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}"

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<style>{CSS}
.sr-only {{ position:absolute; width:1px; height:1px; padding:0; margin:-1px; overflow:hidden; clip:rect(0 0 0 0); white-space:nowrap; border:0; }}
</style>
</head>
<body>
<div class="layout">
  <aside class="sidebar">
    <div class="sb-head">
      <div class="sb-brand">{esc(brand)}</div>
      <div class="sb-title">{esc(subtitle or title)}</div>
      <div class="sb-meta">{meta_html}</div>
    </div>
    <nav class="sb-nav" aria-label="Report sections">{''.join(nav_html)}</nav>
    <div class="sb-foot">{esc(foot)}</div>
  </aside>
  <main class="main">
    <div class="topbar">
      <h1>{esc(title)}</h1>
      {controls}
    </div>
    <div class="content">
      {sections_html}
    </div>
  </main>
</div>
<script>{JS}</script>
</body>
</html>
"""
