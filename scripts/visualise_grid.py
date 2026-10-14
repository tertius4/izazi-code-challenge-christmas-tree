#!/usr/bin/env python3
"""Render every N-tree configuration of a submission as a zoomable HTML grid.

Reuses the evaluator's loading and geometry helpers without modifying it.
"""
import argparse
import csv
import json
import os
import sys
import webbrowser
from decimal import Decimal

import numpy as np
from shapely.ops import unary_union

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from evaluator import (SCALE_FACTOR, calculate_bounding_box_side,  # noqa: E402
                       check_overlap, create_tree_polygon, load_submission)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_SUBMISSION = os.path.join(ROOT, 'data', 'output', 'results.csv')
DEFAULT_OUTPUT = os.path.join(ROOT, 'data', 'output', 'grid.html')
REPORT = os.path.join(ROOT, 'data', 'output', 'report.csv')
SCALE = float(SCALE_FACTOR)


def load_origins(path):
    """Map n -> strategy name from report.csv, if it exists."""
    if not os.path.exists(path):
        return {}
    with open(path, newline='') as f:
        return {int(r['n']): r['origin'] for r in csv.DictReader(f)
                if r['n'].isdigit()}


def build_cell(n, config_df, origin):
    polygons = [create_tree_polygon(r['x'], r['y'], r['deg'])
                for _, r in config_df.iterrows()]
    side = calculate_bounding_box_side(polygons)
    minx, miny, maxx, maxy = (float(Decimal(str(b)) / SCALE_FACTOR)
                              for b in unary_union(polygons).bounds)
    width, height = maxx - minx, maxy - miny
    # Same centred square as evaluator.visualize_configuration.
    sq = float(side)
    sx = minx if width >= height else minx - (sq - width) / 2
    sy = miny if height >= width else miny - (sq - height) / 2
    overlaps, _ = check_overlap(polygons)
    trees = [np.round(np.array(p.exterior.coords) / SCALE, 6).ravel().tolist()
             for p in polygons]
    return {'n': n, 'side': round(sq, 6), 'sq': [round(sx, 6), round(sy, 6)],
            'overlap': bool(overlaps), 'origin': origin, 'trees': trees}


HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Tree Packing Grid</title>
<style>
:root{--bg:#f6f7f9;--fg:#1c2230;--card:#fff;--line:#d5d9e0;--bad:#d62728;--sq:#d62728}
@media (prefers-color-scheme:dark){:root{--bg:#14171d;--fg:#e6e9ef;--card:#1d222b;--line:#333a47}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:14px system-ui,sans-serif;overflow:hidden}
header{position:fixed;inset:0 0 auto 0;height:40px;display:flex;gap:16px;align-items:center;
  padding:0 16px;background:var(--card);border-bottom:1px solid var(--line);z-index:5}
header small{opacity:.7}
#viewport{position:absolute;inset:40px 0 0 0;overflow:hidden;cursor:grab;touch-action:none}
#grid{position:absolute;left:0;top:0;transform-origin:0 0;display:grid;gap:12px;padding:16px;
  grid-template-columns:repeat(10,220px)}
.cell{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:6px;cursor:pointer}
.cell.bad{border-color:var(--bad)}
.cell h4{margin:0 0 4px;font-size:12px;font-weight:600;display:flex;justify-content:space-between}
.cell svg{width:100%;aspect-ratio:1;display:block}
#detail{position:fixed;inset:0;background:rgba(0,0,0,.6);z-index:10;display:none;
  padding:56px 24px 24px}
#detail.open{display:block}
#dbox{position:relative;width:100%;height:100%;background:var(--card);border-radius:8px;
  overflow:hidden;cursor:grab;touch-action:none}
#dsvgwrap{position:absolute;left:0;top:0;transform-origin:0 0}
#dtitle{position:absolute;left:12px;top:10px;z-index:2;font-weight:600;
  background:var(--card);padding:2px 8px;border-radius:4px}
#dhelp{position:fixed;top:20px;left:0;right:0;text-align:center;color:#fff;z-index:11;font-size:12px}
.sq{fill:none;stroke:var(--sq);stroke-dasharray:4 3;vector-effect:non-scaling-stroke;stroke-width:1.5}
.tree{vector-effect:non-scaling-stroke;stroke-width:1;fill-opacity:.5}
</style></head><body>
<header><strong>Tree packing: all configurations</strong>
<small>wheel: zoom &middot; drag: pan &middot; click a cell for detail &middot; double-click: reset</small></header>
<div id="viewport"><div id="grid"></div></div>
<div id="detail"><div id="dhelp">wheel: zoom &middot; drag: pan &middot; &larr;/&rarr;: prev/next N &middot; Esc: close</div>
<div id="dbox"><div id="dtitle"></div><div id="dsvgwrap"></div></div></div>
<script>
const DATA = __DATA__;
const NS = 'http://www.w3.org/2000/svg';
const STOPS = [[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
function color(t){
  const x = t*(STOPS.length-1), i = Math.min(Math.floor(x), STOPS.length-2), f = x-i;
  return 'rgb('+STOPS[i].map((c,k)=>Math.round(c+(STOPS[i+1][k]-c)*f)).join(',')+')';
}
function makeSvg(c){
  const pad = c.side*0.06, svg = document.createElementNS(NS,'svg');
  svg.setAttribute('viewBox', [c.sq[0]-pad, -(c.sq[1]+c.side+pad), c.side+2*pad, c.side+2*pad].join(' '));
  const g = document.createElementNS(NS,'g'); g.setAttribute('transform','scale(1,-1)');
  c.trees.forEach((t,i)=>{
    const p = document.createElementNS(NS,'polygon'), pts = [];
    for(let k=0;k<t.length;k+=2) pts.push(t[k]+','+t[k+1]);
    const col = color(c.n>1 ? i/(c.n-1) : 0);
    p.setAttribute('points', pts.join(' ')); p.setAttribute('class','tree');
    p.setAttribute('fill', col); p.setAttribute('stroke', col); g.appendChild(p);
  });
  const r = document.createElementNS(NS,'rect');
  r.setAttribute('x',c.sq[0]); r.setAttribute('y',c.sq[1]);
  r.setAttribute('width',c.side); r.setAttribute('height',c.side); r.setAttribute('class','sq');
  g.appendChild(r); svg.appendChild(g); return svg;
}
// Generic pan/zoom: applies translate+scale to `target` inside `view`.
function panzoom(view, target, onClick){
  const s = {k:1,x:0,y:0}; let drag = null;
  const apply = () => target.style.transform = `translate(${s.x}px,${s.y}px) scale(${s.k})`;
  view.addEventListener('wheel', e => {
    e.preventDefault();
    const r = view.getBoundingClientRect(), px = e.clientX-r.left, py = e.clientY-r.top;
    const k = Math.min(200, Math.max(0.05, s.k*Math.exp(-e.deltaY*0.0015)));
    s.x = px-(px-s.x)*k/s.k; s.y = py-(py-s.y)*k/s.k; s.k = k; apply();
  }, {passive:false});
  view.addEventListener('pointerdown', e => {
    drag = {x:e.clientX, y:e.clientY, sx:s.x, sy:s.y, moved:false}; view.setPointerCapture(e.pointerId);
    view.style.cursor = 'grabbing';
  });
  view.addEventListener('pointermove', e => {
    if(!drag) return;
    const dx = e.clientX-drag.x, dy = e.clientY-drag.y;
    if(Math.abs(dx)+Math.abs(dy) > 4) drag.moved = true;
    s.x = drag.sx+dx; s.y = drag.sy+dy; apply();
  });
  view.addEventListener('pointerup', e => {
    const d = drag; drag = null; view.style.cursor = 'grab';
    if(d && !d.moved && onClick) onClick(e);
  });
  const reset = fit => { Object.assign(s, fit || {k:1,x:0,y:0}); apply(); };
  view.addEventListener('dblclick', () => reset());
  return {reset};
}
const grid = document.getElementById('grid');
DATA.forEach(c => {
  const d = document.createElement('div'); d.className = 'cell'+(c.overlap?' bad':''); d.dataset.n = c.n;
  d.title = c.origin ? 'origin: '+c.origin : '';
  d.innerHTML = '<h4><span>N = '+c.n+'</span><span>'+c.side.toFixed(4)+'</span></h4>';
  d.appendChild(makeSvg(c)); grid.appendChild(d);
});
const viewport = document.getElementById('viewport');
const gridPZ = panzoom(viewport, grid, e => {
  const cell = e.target.closest('.cell'); if(cell) openDetail(+cell.dataset.n);
});
// Detail overlay
const detail = document.getElementById('detail'), dbox = document.getElementById('dbox');
const wrap = document.getElementById('dsvgwrap'), title = document.getElementById('dtitle');
const detailPZ = panzoom(dbox, wrap); let cur = null;
function openDetail(n){
  cur = n; const c = DATA.find(x => x.n === n);
  wrap.innerHTML = ''; const svg = makeSvg(c);
  const size = Math.min(dbox.clientWidth, dbox.clientHeight);
  svg.setAttribute('width', size); svg.setAttribute('height', size); wrap.appendChild(svg);
  detailPZ.reset({k:1, x:(dbox.clientWidth-size)/2, y:(dbox.clientHeight-size)/2});
  title.textContent = 'N = '+n+'  side = '+c.side.toFixed(6)+(c.overlap?'  OVERLAP':'')+(c.origin?'  ('+c.origin+')':'');
  detail.classList.add('open');
}
detail.addEventListener('click', e => { if(e.target === detail) detail.classList.remove('open'); });
addEventListener('keydown', e => {
  if(!detail.classList.contains('open')) return;
  if(e.key === 'Escape') detail.classList.remove('open');
  else if(e.key === 'ArrowRight' && cur < DATA[DATA.length-1].n) openDetail(cur+1);
  else if(e.key === 'ArrowLeft' && cur > DATA[0].n) openDetail(cur-1);
});
</script></body></html>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('submission', nargs='?', default=DEFAULT_SUBMISSION)
    ap.add_argument('-o', '--output', default=DEFAULT_OUTPUT)
    ap.add_argument('--no-open', action='store_true', help="don't open the browser")
    args = ap.parse_args()

    df = load_submission(args.submission)
    origins = load_origins(REPORT)
    cells = [build_cell(int(n), g, origins.get(int(n), ''))
             for n, g in df.groupby('n')]

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, 'w') as f:
        f.write(HTML.replace('__DATA__', json.dumps(cells, separators=(',', ':'))))
    bad = [c['n'] for c in cells if c['overlap']]
    print(f"Wrote {len(cells)} configurations to {args.output}"
          + (f" (overlaps in N={bad})" if bad else ""))
    if not args.no_open:
        webbrowser.open('file://' + os.path.abspath(args.output))


if __name__ == '__main__':
    main()
