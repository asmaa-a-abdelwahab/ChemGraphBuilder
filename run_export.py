import glob
from schema_relationship_exporter import SchemaExportConfig, SchemaRelationshipExporter

cfg = SchemaExportConfig(
    all_data_connected="Data/AllDataConnected.csv",
    compound_gene_relationship="Data/Relationships/Compound_Gene_Relationship/Compound_Gene_Relationship.csv",
    linkdb_compound_compound_files=glob.glob("Data/Relationships/Cpd_Cpd_CoOccurrence/CID_*.csv"),
    linkdb_gene_symbol_files=glob.glob("Data/Relationships/Cpd_Gene_CoOccurrence/Cpd_Gene_CoOccurrence*.csv"),  # extend to all files
    sdq_consolidatedcompoundtarget_files=glob.glob("Data/Relationships/Compound_Gene_Relationship/Compound_Gene_Interaction_Outside_PubChem*.csv"),  # extend to all files
    compound_transformation="Data/Relationships/Compound_Transformation.csv",
    out_dir="Data/Relationships/SCHEMA",
)

SchemaRelationshipExporter(cfg).export_all()
print("Done. Outputs in Data/Relationships/SCHEMA/")
