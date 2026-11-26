from pathlib import Path
from copy import deepcopy
import xml.etree.ElementTree as E
import json
import cairosvg
NS = 'http://www.w3.org/2000/svg'
E.register_namespace('', NS)
Q = lambda s: '{' + NS + '}' + s
import argparse
ap = argparse.ArgumentParser(description='Reflow compatible external SVG without changing tree paths')
ap.add_argument('source_svg', type=Path)
ap.add_argument('label_json', type=Path)
ap.add_argument('output', type=Path)
args = ap.parse_args()
labels = json.loads(args.label_json.read_text(encoding='utf-8-sig'))
p = args.source_svg
args.output.mkdir(exist_ok=True, parents=True)
old = E.parse(p).getroot()
W, H = (175 * 72 / 25.4, 215 * 72 / 25.4)
r = E.Element(Q('svg'), {'width': '175mm', 'height': '215mm', 'viewBox': f'0 0 {W} {H}'})
E.SubElement(r, Q('rect'), {'width': str(W), 'height': str(H), 'fill': 'white'})
for i in (1, 2):
    r.append(deepcopy(old[i]))

def text(x, y, s, size=7.5, bold=False, **attrs):
    a = {'x': str(x), 'y': str(y), 'font-family': 'Arial', 'font-size': str(size), 'fill': 'black', **attrs}
    if bold:
        a['font-weight'] = 'bold'
    e = E.SubElement(r, Q('text'), a)
    e.text = s
    return e
text(8, 15, 'a', 11, True)
text(23, 15, labels['genome_title'], 9, True)
text(8, 31, labels['host_header'])
text(8, 42, labels['host_span'])
for i in (6, 7):
    e = deepcopy(old[i])
    e.set('y1', '52')
    e.set('y2', '52')
    r.append(e)
text(8, 67, labels['query_span'])
g = E.SubElement(r, Q('g'), {'transform': 'translate(0 12)'})
for i in range(10, 53):
    g.append(deepcopy(old[i]))
text(8, 101, labels['marker_span'])
text(8, 113, labels['marker_note'])
text(266, 15, 'b', 11, True)
text(281, 15, labels['marker_title'], 9, True)
for k in range(6):
    x = 268 + k % 3 * 75
    y = 33 + k // 3 * 31
    text(x, y, ''.join(old[57 + 2 * k].itertext()), 8, True)
    text(x, y + 12, ''.join(old[58 + 2 * k].itertext()), 8)
text(268, 96, labels['marker_total'])
text(268, 107, labels['marker_sources'])
text(268, 118, 'SH-aLRT / UFBoot (%)')
text(8, 137, 'c', 11, True)
text(23, 137, labels['tree_title'], 9, True)
for i in range(123, 131):
    e = deepcopy(old[i])
    if e.tag == Q('text'):
        e.set('y', '153')
        e.set('font-size', '7.5')
        r.append(e)
    else:
        g = E.SubElement(r, Q('g'), {'transform': 'translate(0 24)'})
        g.append(e)
text(247, 153, 'SH-aLRT / UFBoot')
E.SubElement(r, Q('polygon'), {'points': '12,164 19,161 19,167', 'fill': '#D8DDE2', 'stroke': '#9BA4AE', 'stroke-width': '.6'})
text(25, 167, 'Collapsed reference clade')
E.SubElement(r, Q('rect'), {'x': '177', 'y': '161', 'width': '6', 'height': '6', 'fill': 'white', 'stroke': '#9BA4AE', 'stroke-width': '.7'})
text(190, 167, 'Omitted sisters at the same node')
tree = deepcopy(old[75])
tree.set('transform', 'translate(8 152) scale(0.37 0.5625)')
r.append(tree)
for i in range(76, 123):
    e = deepcopy(old[i])
    e.set('y', str(152 + (float(e.attrib['y']) - 112) * 1.25))
    e.set('font-size', '7.5')
    r.append(e)
text(102.7, 569, 'Branch length (amino-acid substitutions/site)', **{'text-anchor': 'middle'})
original_paths = [e.attrib.get('d') for e in old[75].iter(Q('path'))]
assert original_paths == [e.attrib.get('d') for e in tree.iter(Q('path'))]
svg=args.output/'Figure5.svg'; E.ElementTree(r).write(svg,encoding='utf-8',xml_declaration=True)
cairosvg.svg2pdf(url=str(svg),write_to=str(svg.with_suffix('.pdf')))
cairosvg.svg2png(url=str(svg),write_to=str(svg.with_suffix('.png')),dpi=600)
