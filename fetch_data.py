"""
fetch_data.py — Nelson Lab Demo (Alexandra Nelson, UCSF)
Downloads Allen Mouse Brain Atlas ISH expression data for key striatal genes
across motor-circuit brain regions.

Key question: across the motor circuit (dorsal striatum, accumbens, motor cortex,
thalamus), which regions show the co-expression profile that predicts the
dMSN hyperactivation described in Ryan et al. (2024, Cell Reports)?

The Nelson lab (Ryan et al. 2024) found that levodopa-induced dyskinesia (LID)
is driven by a subpopulation of dorsal striatal dMSNs with higher Drd1 expression,
higher dopamine-dependent excitability, and more glutamatergic input from motor cortex.
"""

import urllib.request, json, csv, os

OUT = os.path.dirname(os.path.abspath(__file__))

STRUCTS = {
    "CP":   247,   # Caudoputamen (dorsal striatum — primary LID site)
    "ACB":  56,    # Nucleus accumbens (ventral striatum — reward/motivation)
    "OT":   131,   # Olfactory tubercle (ventral striatum extension)
    "MOp":  985,   # Primary motor cortex (sends glutamatergic drive to dMSNs)
    "Thal": 549,   # Thalamus (another excitatory input to striatum)
    "SNr":  381,   # Substantia nigra reticulata (dMSN output target)
    "Cereb": 512,  # Cerebellum (control region, minimal dopamine/LID involvement)
    "Hipp": 382,   # Hippocampus (comparison)
}

EXPERIMENTS = {
    "Drd1":  80266582,   # D1 dopamine receptor — direct pathway MSN marker
    "Drd2":  77674878,   # D2 dopamine receptor — indirect pathway MSN marker
    "Pdyn":  75929626,   # Prodynorphin — D1-MSN peptide co-transmitter
    "Grin2b": 69257725,  # GluN2B NMDA receptor subunit — excitatory drive indicator
}

print("Querying Allen Mouse Brain Atlas ISH API ...")
results = {name: {} for name in STRUCTS}
for gene, exp_id in EXPERIMENTS.items():
    print(f"  {gene} (exp {exp_id}) ...")
    for name, struct_id in STRUCTS.items():
        url = (f"https://api.brain-map.org/api/v2/data/query.json?"
               f"criteria=model::StructureUnionize,"
               f"rma::criteria,section_data_set[id$eq{exp_id}],"
               f"structure[id$eq{struct_id}]&num_rows=1")
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                d = json.loads(r.read())
            msg = d.get("msg", [])
            val = float(msg[0].get("expression_density", 0) or 0) if msg else 0.0
            results[name][gene] = round(val, 6)
        except Exception as e:
            print(f"    WARNING {name}: {e}")
            results[name][gene] = 0.0

with open(os.path.join(OUT, "striatum_genes.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["region"] + list(EXPERIMENTS.keys()))
    for name in STRUCTS:
        w.writerow([name] + [results[name][g] for g in EXPERIMENTS])

print("\nWrote striatum_genes.tsv")
print("\nKey values (Drd1/Drd2/Grin2b):")
for name in STRUCTS:
    r = results[name]
    ratio = r['Drd1'] / max(r['Drd2'], 1e-6)
    print(f"  {name:<6}: Drd1={r['Drd1']:.4f}  Drd2={r['Drd2']:.4f}  "
          f"D1/D2={ratio:.2f}  Grin2b={r['Grin2b']:.4f}")
