#!/usr/bin/env python3
import argparse, csv, html, statistics
from pathlib import Path

def rows(path):
    with Path(path).open() as h: return list(csv.DictReader(h, delimiter="\t"))

def paf(path, start, end, min_aligned=0, min_identity=.90, min_mapq=20):
    """Every alignment passing the filters, not one per contig.

    A contig that aligns to two separate island segments contributes both: the
    island is repeat-rich and targeted assemblies are fragmented, so collapsing
    to a single best alignment per contig discards real reference coverage.
    Use best_per_contig() where a per-contig view is wanted.

    T1: min_aligned defaults to 0 because coverage must not have a length floor.
    A 213 bp exact match is 213 bp of island recovered whether or not the contig
    it came from counts as "supporting evidence". Identity and MAPQ do the
    specificity work. The floor belongs to supporting_contigs(), below.
    """
    hits=[]
    for line in Path(path).read_text().splitlines():
        x=line.split("\t")
        if len(x)<12: continue
        alen=int(x[10]); ident=int(x[9])/alen if alen else 0; ts=int(x[7]); te=int(x[8]); mapq=int(x[11])
        if te<=start or ts>=end or alen<min_aligned or ident<min_identity or mapq<min_mapq: continue
        hits.append(dict(contig=x[0],qlen=int(x[1]),strand=x[4],start=ts,end=te,aligned=alen,identity=ident,mapq=mapq))
    return sorted(hits, key=lambda x:(x["start"],-x["aligned"]))

def best_per_contig(hits):
    """One representative alignment per contig, for counting supporting contigs."""
    best={}
    for hit in hits:
        old=best.get(hit["contig"])
        if old is None or (hit["aligned"],hit["identity"],hit["mapq"])>(old["aligned"],old["identity"],old["mapq"]): best[hit["contig"]]=hit
    return sorted(best.values(), key=lambda x:(x["start"],-x["aligned"]))

def supporting_contigs(hits, island_start, island_end, min_total_aligned):
    """Contigs carrying at least min_total_aligned bases of ISLAND overlap, unioned
    per contig so an alignment block spanning the island boundary, or two blocks
    from the same contig overlapping each other, isn't double-counted.

    T1: summed per contig, not per alignment block. The island is repeat-rich and
    targeted assemblies are fragmented, so one contig routinely matches in several
    pieces; judging each block on its own discarded contigs that carry plenty of
    island between them. clbR is 213 bp and clbE 249 bp, so a floor above ~200
    can also reject a contig that covers a whole gene.

    Ludmil, revised report finding 3: this used to sum each hit's full aligned-block
    length within the wider +/-10kb recruitment window `paf()` is called over, not
    the portion actually overlapping the canonical island. A contig aligning mostly
    or entirely to a flank -- with little or no real island overlap -- could still
    clear min_total_aligned and be reported as island-supporting. union_length()
    was always clipped to the island correctly; only this per-contig count was not,
    so it now reuses the same clipping, just scoped to one contig's hits at a time.
    """
    by_contig = {}
    for hit in hits:
        by_contig.setdefault(hit["contig"], []).append(hit)
    return sorted(
        contig for contig, contig_hits in by_contig.items()
        if union_length(contig_hits, island_start, island_end) >= min_total_aligned
    )


def union_length(hits, start, end):
    spans=sorted((max(start,h["start"]),min(end,h["end"])) for h in hits if h["end"]>start and h["start"]<end)
    total=0; right=start
    for left,stop in spans:
        total += max(0,stop-max(left,right)); right=max(right,stop)
    return total

def depths(path, contig, start, end):
    out=[]
    with Path(path).open() as h:
        for line in h:
            x=line.rstrip().split("\t")
            if len(x)!=3 or x[0]!=contig: continue
            pos=int(x[1])-1
            if pos<start: continue
            if pos>=end: break
            out.append((pos,int(x[2])))
    return out

def genes(path, contig, start, end):
    out=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"): continue
        x=line.split("\t")
        if len(x)<9 or x[0]!=contig or x[2]!="gene": continue
        left,right=int(x[3])-1,int(x[4]); attrs=dict(v.split("=",1) for v in x[8].split(";") if "=" in v)
        if right>start and left<end: out.append((left,right,x[6],attrs.get("Name",attrs.get("gene","gene"))))
    return out

def structural_call(mh, ms):
    lo,hi=min(mh,ms),max(mh,ms)
    if hi>=.90: call="complete_island"
    elif hi>=.75: call="near_complete_island"
    elif hi>=.20: call="fragmented_island"
    elif hi>=.05: call="partial_island"
    elif hi>0: call="limited_contig_support"
    else: call="read_evidence_only"
    if hi==0: agreement="no_contig_support"
    elif lo==0: agreement="single_assembler_only"
    elif abs(mh-ms)<=.20: agreement="concordant"
    else: agreement="discordant"
    return call,agreement

def main():
    p=argparse.ArgumentParser()
    for name in ("sample","megahit-contigs","metaspades-contigs","megahit-paf","metaspades-paf","raw-depth","read-evidence","gff","output-tsv","output-svg"): p.add_argument(f"--{name}",required=True)
    p.add_argument("--contig",default="NC_017628.1"); p.add_argument("--region-start",type=int,default=2183826); p.add_argument("--region-end",type=int,default=2254594)
    p.add_argument("--island-start",type=int,default=2193826); p.add_argument("--island-end",type=int,default=2244594)
    p.add_argument("--min-aligned-bp",type=int,default=200,
                   help="island bases a contig must carry IN TOTAL to be counted as a "
                        "supporting contig. Applies to the contig count only -- coverage, "
                        "and therefore assembler_agreement and the structural call, count "
                        "every alignment passing identity and MAPQ. Default 200 sits below "
                        "the shortest clb gene (clbR, 213 bp), so a contig covering a whole "
                        "gene is never rejected")
    a=p.parse_args(); ev=rows(a.read_evidence)[0]; length=a.island_end-a.island_start
    # T1: no length floor here. Coverage, and every call derived from it, counts
    # every alignment that passes identity and MAPQ.
    hits={"megahit":paf(a.megahit_paf,a.region_start,a.region_end),"metaspades":paf(a.metaspades_paf,a.region_start,a.region_end)}
    best={k:best_per_contig(v) for k,v in hits.items()}
    supporting={k:supporting_contigs(v,a.island_start,a.island_end,a.min_aligned_bp) for k,v in hits.items()}
    bp={k:union_length(v,a.island_start,a.island_end) for k,v in hits.items()}; cov={k:v/length for k,v in bp.items()}
    call,agreement=structural_call(cov["megahit"],cov["metaspades"])
    fields=["sample","read_evidence","pks_reads","clb_genes_detected","island_breadth_1x","island_breadth_2x","island_breadth_3x","megahit_reference_covered_bp","megahit_reference_coverage","megahit_supporting_contigs","metaspades_reference_covered_bp","metaspades_reference_coverage","metaspades_supporting_contigs","assembler_agreement","final_structural_evidence"]
    result=dict(ev,megahit_reference_covered_bp=bp["megahit"],megahit_reference_coverage=f'{cov["megahit"]:.6f}',megahit_supporting_contigs=len(supporting["megahit"]),metaspades_reference_covered_bp=bp["metaspades"],metaspades_reference_coverage=f'{cov["metaspades"]:.6f}',metaspades_supporting_contigs=len(supporting["metaspades"]),assembler_agreement=agreement,final_structural_evidence=call)
    with Path(a.output_tsv).open("w",newline="") as h: w=csv.DictWriter(h,fields,delimiter="\t");w.writeheader();w.writerow(result)
    W,H,L,R=1600,1100,100,35; PW=W-L-R; sx=lambda x:L+(x-a.region_start)/(a.region_end-a.region_start)*PW
    svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><rect width="100%" height="100%" fill="white"/>','<style>text{font-family:Arial,sans-serif;fill:#202124}.title{font-size:24px;font-weight:bold}.label{font-size:14px;font-weight:bold}.small{font-size:10px}.sub{font-size:13px}</style>',f'<text x="{L}" y="34" class="title">Canonical PKS-island contig validation</text>',f'<text x="{L}" y="58" class="sub">{html.escape(a.sample)} • IHE3034 {a.contig} • {ev["read_evidence"]} read evidence • final: {call}</text>',f'<text x="{L}" y="88" class="label">A  Standard reference clb genes</text>',f'<line x1="{L}" y1="112" x2="{W-R}" y2="112" stroke="#444"/>']
    for left,right,strand,name in genes(a.gff,a.contig,a.region_start,a.region_end):
        x1,x2=sx(left),sx(right); col="#3A86FF" if name.startswith("clb") else "#ADB5BD"; svg.append(f'<rect x="{x1:.1f}" y="101" width="{max(1,x2-x1):.1f}" height="22" fill="{col}" stroke="white"/>')
        if name.startswith("clb"): svg.append(f'<text x="{(x1+x2)/2:.1f}" y="140" font-size="9" transform="rotate(48 {(x1+x2)/2:.1f} 140)">{html.escape(name)}</text>')
    y=180
    for assembler,color,label in (("megahit","#277DA1","B  MEGAHIT"),("metaspades","#F8961E","C  metaSPAdes")):
        svg.append(f'<text x="{L}" y="{y}" class="label">{label}: {len(supporting[assembler])} supporting contigs (&#8805;{a.min_aligned_bp} bp each); {100*cov[assembler]:.1f}% island coverage from {len(best[assembler])} aligned contigs</text>'); y+=18
        if not best[assembler]: svg.append(f'<text x="{L+20}" y="{y+12}" class="sub">No contig aligns at ≥90% identity, MAPQ ≥20</text>'); y+=35
        else:
            for hit in best[assembler][:12]:
                x1,x2=sx(max(a.region_start,hit["start"])),sx(min(a.region_end,hit["end"])); svg.append(f'<rect x="{x1:.1f}" y="{y}" width="{max(1,x2-x1):.1f}" height="13" fill="{color}" stroke="#17212b"/>'); svg.append(f'<text x="{min(W-330,x2+5):.1f}" y="{y+10}" class="small">{html.escape(hit["contig"])} • {hit["aligned"]} bp • {100*hit["identity"]:.1f}%</text>'); y+=18
        y+=12
    dep=depths(a.raw_depth,a.contig,a.region_start,a.region_end); top=max(y+25,500); bottom=top+115; cap=max((d for _,d in dep),default=1)
    svg.append(f'<text x="{L}" y="{top-12}" class="label">D  Raw reference depth (MAPQ ≥40; blank = 0)</text>')
    bins={}
    for pos,d in dep: bins.setdefault(int((pos-a.region_start)/(a.region_end-a.region_start)*PW),[]).append(d)
    for b,values in bins.items():
        height=statistics.mean(values)/cap*115; svg.append(f'<rect x="{L+b}" y="{bottom-height:.1f}" width="1" height="{height:.1f}" fill="#F94144"/>')
    svg += [f'<line x1="{L}" y1="{bottom}" x2="{W-R}" y2="{bottom}" stroke="#444"/>',f'<rect x="{L}" y="{bottom+18}" width="{PW}" height="55" rx="7" fill="#F8F9FA" stroke="#DADCE0"/>',f'<text x="{L+14}" y="{bottom+40}" class="label">{call} • {agreement}</text>',f'<text x="{L+14}" y="{bottom+60}" class="sub">Breadth ≥1×/≥2×/≥3×: {100*float(ev["island_breadth_1x"]):.1f}% / {100*float(ev["island_breadth_2x"]):.1f}% / {100*float(ev["island_breadth_3x"]):.1f}%</text>','</svg>']
    Path(a.output_svg).write_text("\n".join(svg))

if __name__=="__main__": main()
