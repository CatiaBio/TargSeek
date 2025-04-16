# QuickGO Protein Annotation Data

This folder contains protein annotation data downloaded from [QuickGO](https://www.ebi.ac.uk/QuickGO/annotations). Specific Gene Ontology (GO) terms were selected to identify and retrieve proteins located within defined bacterial cell compartments, including:

- Cell outer membrane
- Peptidoglycan-based cell wall
- Gram-negative-bacterium-type cell wall
- External encapsulating structures
- Structural components associated with the bacterial cell wall

## Purpose

The dataset is intended to facilitate bioinformatics analyses of bacterial proteins localized in specific cellular structures, potentially identifying suitable targets for antibody development or other biomedical applications.

## Data Source

Data was retrieved from the European Bioinformatics Institute (EBI) QuickGO resource:

- [https://www.ebi.ac.uk/QuickGO/annotations](https://www.ebi.ac.uk/QuickGO/annotations)

## Data Selection Criteria

Proteins were selected based on Gene Ontology (GO) terms relevant to bacterial cellular compartments and structures, using the following GO relationships: "is_a", "part_of", and "occurs_in". The selected GO terms include:
GO:0046930    pore complex
GO:0015288    porin activity
GO:0090729    toxin activity
GO:0016020    membrane
GO:0005576    extracellular region
GO:0007155    cell adhesion
GO:0030288 – Outer membrane-bounded periplasmic space
GO:0050840 - extracellular matrix binding
GO:0005576 - extracellular region
GO:0005615 - extracellular space
GO:0005618 – Cell wall: A rigid structure outside the plasma membrane, primarily composed of peptidoglycan in bacteria.​
GO:0005886 – Plasma membrane: Semipermeable membrane enclosing the cytoplasm of bacterial cells.​
GO:0009270 – Bacterial-type flagellum motor: Rotary motor driving the rotation of the bacterial flagellum.​
GO:0009271 – Bacterial-type flagellum filament: Long, helical filament that propels the bacterium.​
GO:0009272 – Bacterial-type flagellum hook: Flexible coupling between the basal body and the filament of the flagellum.​
GO:0009273 – Gram-positive cell wall: Thick peptidoglycan-rich wall characteristic of Gram-positive bacteria.​
GO:0009274 – Peptidoglycan-based cell wall: A structural component unique to bacteria, providing rigidity and shape.​
GO:0009275 – Gram-negative cell wall: Thin peptidoglycan layer located between the inner and outer membranes in Gram-negative bacteria.​
GO:0009277 – Bacterial-type flagellum basal body: Anchoring structure of the flagellum embedded in the cell envelope.​
GO:0009278 – Bacterial-type flagellum: Long, whip-like structures enabling bacterial motility.​
GO:0009279 – Cell outer membrane: The outermost layer in Gram-negative bacteria, containing lipopolysaccharides.​
GO:0009288 – Fimbriae: Thin, filamentous structures aiding in bacterial attachment to surfaces.​
GO:0009289 – Pilus: Hair-like appendages on bacterial surfaces involved in adhesion and conjugation.​
GO:0009425 – Lipopolysaccharide biosynthetic process: Formation of lipopolysaccharides, key components of the outer membrane in Gram-negative bacteria.​
GO:0010339 – Cell wall polysaccharide biosynthetic process: Synthesis of polysaccharides that are integral to the bacterial cell wall.​
GO:0019861 – Outer membrane: External membrane found in Gram-negative bacteria, distinct from the plasma membrane.​
GO:0030312 – External encapsulating structure: Structures like capsules or sheaths that enclose the cell, external to the cell wall.​
GO:0030313 – Cell envelope: Combined structure of the plasma membrane, cell wall, and, in Gram-negative bacteria, the outer membrane.​
GO:0097313 – Bacterial-type cell wall organization: Processes involved in the assembly and arrangement of the bacterial cell wall.​
GO:0099008 – Bacterial-type flagellum-dependent cell motility: Movement of bacterial cells facilitated by the rotation of flagella.

## Data Format

Annotations from QuickGO are provided in a tab-delimited (TSV) or CSV format with the following key columns:

- GENE PRODUCT DB  
- GENE PRODUCT ID  
- SYMBOL  
- QUALIFIER  
- GO TERM  
- GO NAME  
- ECO ID  
- GO EVIDENCE CODE  
- REFERENCE  
- WITH/FROM  
- TAXON ID  
- ASSIGNED BY  
- ANNOTATION EXTENSION  
- GO ASPECT

## Usage

This dataset is intended exclusively for research purposes. Please cite QuickGO in any publications or analyses that utilize this data:

> UniProt Consortium. QuickGO: a web-based tool for Gene Ontology searching. [https://www.ebi.ac.uk/QuickGO](https://www.ebi.ac.uk/QuickGO)

## Date of Retrieval
2025-04-10

