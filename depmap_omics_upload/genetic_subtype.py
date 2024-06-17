import pandas as pd
from taigapy import TaigaClient
tc = TaigaClient()


def subtype_egfr_exon19del(maf, subtype_table, colname="EGFR exon 19 del", gene_colname = "HugoSymbol", protein_change_colname = "ProteinChange"):
    exon19_start = 55174722
    exon19_end = 55174820
    subtype_table.loc[maf[(maf[gene_colname] == "EGFR") & 
                          (maf.Pos > exon19_start) & 
                          (maf.Pos < exon19_end) & 
                          (maf[protein_change_colname].str.contains("del"))].ModelID.unique().tolist(), 
                      colname] = True
    return subtype_table

def generate_subtype_matrix(maf, models, single_protein_changes, single_position, hotspot_genes, damaging_genes, gene_colname = "HugoSymbol", protein_change_colname = "ProteinChange", hotspot_colname = "Hotspot", damaging_colname = "LikelyLoF"):
    # initialize matrix
    egfr_colname = "EGFR exon 19 del"
    colnames = ([p[0] + " " + p[1] for p in single_protein_changes] + 
                [p[0] + " " + p[1] + "*" for p in single_position] + 
                [g + " Hotspot" for g in hotspot_genes] + 
                [g + " Damaging" for g in damaging_genes] + 
                [egfr_colname])
    subtype_table = pd.DataFrame(index=models["ModelID"].unique(), columns=colnames)
    subtype_table.loc[set(maf["ModelID"]), colnames] = False
    
    for (gene, pc) in single_protein_changes:
        subtype_table.loc[maf[(maf[gene_colname] == gene) & (maf[protein_change_colname] == pc)].ModelID.tolist(), gene + " " + pc] = True
        
    for (gene, pos) in single_position:
        subtype_table.loc[maf[(maf[gene_colname] == gene) & (maf[protein_change_colname].str.startswith(pos))].ModelID.unique().tolist(), gene + " " + pos + "*"] = True
    
    # add hotspot and damaging
    for gene in hotspot_genes:
        subtype_table.loc[maf[(maf[gene_colname] == gene) & (maf[hotspot_colname] == True)].ModelID.tolist(), gene + " Hotspot"] = True
    
    for gene in damaging_genes:
        subtype_table.loc[maf[(maf[gene_colname] == gene) & (maf[damaging_colname] == True)].ModelID.tolist(), gene + " Damaging"] = True
    
    
    # add EGFR column with more complex logic
    subtype_table = subtype_egfr_exon19del(maf, subtype_table, colname=egfr_colname, gene_colname=gene_colname, protein_change_colname=protein_change_colname)
    
    return subtype_table



def __main__():

    taiga_virtual = {"internal": "internal-24q2-3719",
                    "dmc": "dmc-24q2-5194",
                    "public": "public-24q2-356f"}

    upload_dataset = "genetic-subtypes-7e8a"
    mutation_table_fn = "OmicsSomaticMutations"
    model_metadata_fn = "Model"

    # list of subtypes defined by one gene and one protein change
    single_protein_changes = [("KRAS", "p.G12D"), ("KRAS", "p.G12C"),
                            ("BRAF", "p.V600E"), ("EGFR", "p.L858R"), ("JAK2", "V617F")]

    # list of subtypes defined by one gene and one position, including all possible protein changes
    single_position = [("KRAS", "p.G12"), ("KRAS", "p.G13"), ("KRAS", "p.Q61"),
                    ("NRAS", "p.G12"), ("NRAS", "p.G13"), ("NRAS", "p.Q61"),
                    ("HRAS", "p.G12"), ("HRAS", "p.G13"), ("HRAS", "p.Q61"),
                    ("PIK3CA", "p.E542"), ("PIK3CA", "p.E545"), ("PIK3CA", "p.H1047"),
                    ]

    # list of subtypes defined by gene and hotspot == True
    hotspot_genes = ["TP53", "NF1"]

    # list of subtypes defined by gene and damaging == True
    damaging_genes = ["TP53", "NF1"]

    maf = tc.get(name=taiga_virtual['internal'], file=mutation_table_fn)
    model = tc.get(name=taiga_virtual['internal'], file=model_metadata_fn)

    t = generate_subtype_matrix(maf, model, single_protein_changes, single_position, hotspot_genes, damaging_genes)

    t.to_csv("subtypes.csv")

    tc.update_dataset(
        dataset_id=upload_dataset,
        changes_description="added hotspot and damaging",
        upload_files=[
            {
                "path": "subtypes.csv",
                "name": "genetic_subtypes",
                "format": "TableCSV",
                "encoding": "utf-8",
            },
        ],
        add_all_existing_files=True,
    )