"""MathML definitions for the manuscript's editable equations."""
import json
from pathlib import Path
from xml.sax.saxutils import escape


def el(tag, *parts):
    return '<' + tag + '>' + ''.join(parts) + '</' + tag + '>'


def mi(value): return el('mi', escape(str(value)))
def label(value): return '<mi mathvariant="normal">' + escape(str(value)) + '</mi>'
def mn(value): return el('mn', str(value))
def mo(value): return el('mo', escape(value))
def row(*parts): return el('mrow', *parts)
def sub(base, index): return el('msub', base, index)
def sup(base, power): return el('msup', base, power)
def frac(top, bottom): return el('mfrac', top, bottom)
def brackets(value, left='(', right=')'): return row(mo(left), value, mo(right))
def text(value): return el('mtext', escape(value))
def sum_under(condition, expression): return row(el('munder', mo('∑'), condition), expression)
def absval(value): return brackets(value, '|', '|')
def indicator(condition): return row(mi('𝟙'), brackets(condition, '[', ']'))
def comma(): return row(mo(','), text('  '))
def op(value): return row(text(' '),mo(value),text(' '))
def index(value): return row(*(mi(c) if c.isalpha() else mo(c) for c in value))
def math(value): return '<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">'+value+'</math>'


def definitions():
    specs = {}
    def add(key, number, expression): specs[key]={'number':str(number),'mathml':math(expression)}
    add('M1',1,row(sub(mi('I'),index('uv')),op('='),frac(sub(mi('M'),index('uv')),sub(mi('A'),index('uv')))))

    directed=lambda symbol, direction:sub(mi(symbol),index(direction))
    boundary=lambda side:sub(mi('B'),row(label(side),mo(','),index('u→v')))
    first=row(directed('D','u→v'),op('='),directed('R','uv'),boundary('L'),boundary('R'))
    second=row(sub(mi('S'),brackets(index('u,v'),'{','}')),op('='),mi('max'),brackets(row(directed('D','u→v'),mo(','),directed('D','v→u'))))
    add('M2',2,row(first,comma(),second))

    edge=mi('e'); pk=lambda a: row(sub(mi('p'),mi(a)),brackets(edge))
    eset=lambda a: sup(sub(mi('E'),mi(a)),mi('k'))
    both=row(edge,mo('∈'),eset('A'),mo('∩'),eset('B'))
    numerator=row(mn(2),sum_under(both,row(mi('min'),brackets(row(pk('A'),mo(','),pk('B'))))))
    denominator=row(sum_under(row(edge,mo('∈'),eset('A')),pk('A')),mo('+'),sum_under(row(edge,mo('∈'),eset('B')),pk('B')))
    add('M3',3,row(sub(mi('O'),mi('k')),op('='),frac(numerator,denominator)))

    density=lambda idx:sub(mi('d'),index(idx))
    count=lambda idx:sub(mi('C'),index(idx))
    share=lambda idx:sub(mi('s'),index(idx))
    length=sub(mi('L'),mi('u'))
    sumdensity=row(el('munderover',mo('∑'),row(mi('v'),mo('='),mn(1)),mn(4)),density('vct'))
    add('M4',4,row(density('uct'),op('='),frac(count('uct'),length),comma(),share('uct'),op('='),frac(density('uct'),sumdensity)))
    delta=sub(mi('Δ'),index('uc'))
    add('M5',5,row(delta,brackets(row(sub(mi('t'),mi('a')),mo(','),sub(mi('t'),mi('b')))),op('='),sub(mi('s'),row(mi('u'),mi('c'),sub(mi('t'),mi('b')))),mo('−'),sub(mi('s'),row(mi('u'),mi('c'),sub(mi('t'),mi('a'))))))
    lower=row(el('munder',mi('min'),mi('c')),delta,mo('>'),mn(0))
    upper=row(el('munder',mi('max'),mi('c')),delta,mo('<'),mn(0))
    magnitude=row(el('munder',text('median'),mi('c')),absval(delta),mo('≥'),mn('0.05'))
    add('M6',6,row(sub(mi('H'),mi('u')),op('='),indicator(row(lower,mo('∨'),upper)),mo('·'),indicator(magnitude)))

    rank=lambda symbol:sub(mi(symbol),mi('i'))
    mean=lambda symbol:el('mover',mi(symbol),mo('¯'))
    centered=lambda symbol:brackets(row(rank(symbol),mo('−'),mean(symbol)))
    numerator=sum_under(mi('i'),row(centered('R'),centered('Q')))
    denominator=el('msqrt',row(sum_under(mi('i'),sup(centered('R'),mn(2))),sum_under(mi('i'),sup(centered('Q'),mn(2)))))
    add('S1','S2',row(sub(mi('ρ'),label('s')),op('='),frac(numerator,denominator)))
    rate=lambda a:frac(sub(mi('C'),index(a+'l')),sub(mi('L'),mi(a)))
    add('S2','S3',row(sup(sub(text('TPM'),index('fl')),mo('*')),op('='),sup(mn(10),mn(6)),frac(rate('f'),sum_under(row(mi('g'),mo('∈'),mi('F')),rate('g')))))
    indicator_null=indicator(row(sub(mi('T'),mi('b')),mo('≥'),sub(mi('T'),text('obs'))))
    numerator=row(mn(1),mo('+'),el('munderover',mo('∑'),row(mi('b'),mo('='),mn(1)),mi('B')),indicator_null)
    add('S3','S1',row(sub(mi('p'),text('add-one')),op('='),frac(numerator,row(mi('B'),mo('+'),mn(1))),comma(),mi('B'),op('='),mn(64)))
    return specs


if __name__=='__main__':
    Path(__file__).with_name('equations.json').write_text(json.dumps(definitions(),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
