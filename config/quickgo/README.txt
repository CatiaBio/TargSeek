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

🧫 Cell Surface & Extracellular
GO ID	Name	Description
GO:0005576	Extracellular region	Any space outside the cell, including secreted proteins and fluid-exposed molecules.
GO:0005615	Extracellular space	The space external to the cell membrane, where secreted proteins may be detected.
GO:0005618	Cell wall	Rigid structure outside the plasma membrane; varies between Gram-positive and Gram-negative bacteria.
GO:0019861	Outer membrane	The outermost membrane in Gram-negative bacteria, enriched in LPS and porins.
GO:0030312	External encapsulating structure	Surface structures such as capsules or sheaths that enclose the bacterial cell.
GO:0030313	Cell envelope	The multilayered structure enclosing the cell, including membranes and the wall.
GO:0030288	Outer membrane-bounded periplasmic space	Space between the inner and outer membranes in Gram-negative bacteria.
GO:0016020	Membrane	Generic term for any lipid bilayer surrounding a compartment, including plasma, inner, and outer membranes.

🧬 Motility (Flagella)
GO ID	Name	Description
GO:0009278	Bacterial-type flagellum	The full motility structure (motor, hook, filament).
GO:0009270	Flagellum motor	The rotary engine embedded in the cell envelope driving flagellar rotation.
GO:0009271	Flagellum filament	The long, whip-like filament that extends into the extracellular space.
GO:0009272	Flagellum hook	Connects the motor and filament; provides flexibility and torque transfer.
GO:0009277	Flagellum basal body	The anchoring structure embedded in the membrane layers.
GO:0099008	Flagellum-dependent cell motility	Movement of cells using rotation of flagella.

🦠 Adhesion & Appendages
GO ID	Name	Description
GO:0007155	Cell adhesion	Proteins that mediate attachment to host cells or surfaces.
GO:0050840	Extracellular matrix binding	Proteins that bind host ECM components like collagen or fibronectin.
GO:0009288	Fimbriae	Short, hair-like structures aiding in bacterial adhesion.
GO:0009289	Pilus	Long filamentous structures used for adhesion or conjugation.

🧱 Cell Wall & Membrane Types
GO ID	Name	Description
GO:0009273	Gram-positive cell wall	Thick peptidoglycan layer with teichoic acids; lacks outer membrane.
GO:0009274	Peptidoglycan-based cell wall	Common to all bacteria; provides rigidity and shape.
GO:0009275	Gram-negative cell wall	Thin peptidoglycan layer between inner and outer membranes.
GO:0009276	Gram-negative-bacterium-type cell wall	Specific term for Gram-negative wall structure, similar to above.
GO:0009279	Cell outer membrane	Specific outer membrane found in Gram-negative bacteria.

🧬 Transport & Permeability
GO ID	Name	Description
GO:0015288	Porin activity	Channel proteins forming pores in the outer membrane.
GO:0046930	Pore complex	Multi-protein assemblies that form selective pores across membranes.

💥 Toxins & Biosynthesis
GO ID	Name	Description
GO:0090729	Toxin activity	Proteins that damage or disrupt host functions.
GO:0009425	Lipopolysaccharide biosynthetic process	Assembly of LPS — a major antigen in Gram-negative bacteria.
GO:0010339	Cell wall polysaccharide biosynthesis	Creation of complex sugars for the cell wall and capsule.
GO:0097313	Cell wall organization	Processes that build, maintain, or remodel the bacterial cell wall.

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

