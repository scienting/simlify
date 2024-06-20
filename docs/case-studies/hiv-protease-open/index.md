# HIV protease

This guide focuses on using Amber to simulate HIV-1 protease, a crucial enzyme in the life cycle of the human immunodeficiency virus (HIV) and a prime target for antiretroviral therapy.
The following sections will cover the practical aspects of using Amber for HIV-1 protease simulations, including system preparation, force field selection, simulation protocols, and data analysis techniques.

## Background

HIV-1 protease plays a critical role in the maturation of HIV particles, making it an essential component in the virus's replication process.
HIV-1 protease enables the formation of infectious viral particles by cleaving viral polyproteins into functional proteins.
Understanding this enzyme's structure, dynamics, and function is paramount for developing effective inhibitors and advancing HIV treatment strategies.

The significance of HIV-1 protease extends beyond its immediate role in viral replication:

1.  **Public health impact:** HIV/AIDS remains a global health challenge, affecting millions of people worldwide. Insights gained from studying HIV-1 protease contribute to the ongoing efforts to combat this pandemic.
2.  **Drug design:** As a prime target for antiretroviral drugs, detailed knowledge of HIV-1 protease's structure and dynamics is crucial for rational drug design and developing new therapeutic strategies.
3.  **Resistance mechanisms:** HIV's rapid mutation rate often leads to drug resistance. MD simulations can help elucidate the molecular basis of these resistance mechanisms, informing the design of more robust inhibitors.

From a biophysical perspective, HIV-1 protease presents several exciting characteristics that make it an excellent subject for MD simulations:

1.  **Structural flexibility:** The enzyme exhibits significant conformational changes during its catalytic cycle, particularly in the "flap" regions that control access to the active site.
    MD simulations can capture these dynamic processes, providing insights into the enzyme's mechanism.
2.  **Homodimeric nature:** HIV-1 protease functions as a homodimer, with the active site formed at the interface between two identical subunits.
    This symmetry adds an interesting dimension to the simulation setup and analysis.
3.  **Substrate specificity:** The enzyme recognizes and cleaves specific sequences in viral polyproteins.
    MD simulations can help elucidate the molecular basis of this specificity and how mutations might affect it.
4.  **Water-Mediated Interactions:** Water molecules play a crucial role in the enzyme's catalytic mechanism and mediating protein-ligand interactions.
 Explicit solvent MD simulations are particularly valuable for studying these effects.
5.  **Allosteric effects:** Recent studies have suggested the presence of allosteric sites in HIV-1 protease, opening new avenues for drug design.
    MD simulations can help identify and characterize these sites.

By running Amber simulations on HIV-1 protease, researchers can gain valuable insights into these biophysical characteristics, contributing to fundamental science and applied research in drug discovery.

Understood. I appreciate the flexibility to reorganize and expand on the information. Here's a revised version of the section, incorporating your notes while adding some additional relevant information:

## System preparation

System preparation is a crucial step in molecular dynamics (MD) simulations that significantly impacts the quality and reliability of results.
This process involves selecting an appropriate starting structure and setting up the system to accurately represent physiological conditions.
Proper preparation minimizes artifacts and ensures that simulations reflect the true behavior of the protein in its biological context.

### Protein selection

By meticulously selecting a starting structure based on these criteria, we establish a strong foundation for our MD simulations of HIV-1 protease.
This careful selection process increases the likelihood of obtaining biologically relevant and computationally robust results, which are crucial for understanding the protein's behavior and for applications such as drug design and the study of resistance mechanisms.
Remember that the chosen structure will serve as the basis for all subsequent steps in the simulation process, including system setup, energy minimization, and production runs.
Therefore, the time invested in selecting an appropriate structure is well spent and can save considerable effort and computational resources in the long run.

We will use the Protein Data Bank (PDB) as our primary source of structural data.
When selecting a structure, consider the following criteria:

1.  **Resolution:** Prioritize high-resolution structures (typically < 2.0 Å) for more accurate atomic positions.
    Lower resolution (i.e., > 2.0 Å) structures may lead to simulation inaccuracies due to less precise atomic coordinates.
2.  **Validation Scores:** Use the PDB validation reports to assess structure quality.
    Key metrics include:
    -   Clashscore: Lower values indicate fewer steric clashes between atoms.
    -   Ramachandran outliers: Fewer outliers suggest better backbone geometry.
    -   Rotamer outliers: Fewer outliers indicate more reliable side-chain conformations.
    -   Overall quality at a glance: Provides a quick assessment of the structure's quality.
3.  **Completeness:** Choose structures with minimal missing atoms or residues to reduce uncertainties in the simulation.
4.  **Apo structure:** Use a ligand-free (apo) structure to study the intrinsic dynamics of HIV-1 protease without bias from bound ligands or inhibitors.
5.  **HIV variant:** Select a structure representing the most common HIV variant, such as HIV-1 subtype B, which is prevalent in North America, Western Europe, and Australia.
6.  **Experimental method:** While X-ray crystallography is common, consider high-quality structures from other methods like cryo-electron microscopy (cryo-EM) or NMR spectroscopy if they offer advantages in physiological relevance or completeness.
7.  **Publication date and citations:** Balance recent structures benefiting from improved techniques with well-established, highly cited structures that allow comparison with previous studies.
8.  **Physiological relevance:** Prefer structures determined under conditions mimicking physiological environments, such as appropriate pH and temperature.
9.  **Authors and laboratory reputation:** Structures from laboratories with expertise in HIV protease can provide additional confidence in data quality.

To begin the search, we used the following parameters in the PDB:

-   **Full Text**: `HIV Protease`
-   **Polymer Entity Type** is `Protein`
-   **Refinement Resolution** is between `0.5` to `2` (upper included).
-   **Enzyme Classification Name** is `Hydrolases`.
-   **Scientific Name of the Source Organism** is `Human immunodeficiency virus 1`.

!!! note

    It's worth noting that while we aim for an apo structure to avoid bias, in some cases, the highest quality available structures may be ligand-bound.
    If this is the case, we'll need to carefully remove the ligand and consider running a short equilibration simulation to allow the protein to relax into its unbound state.

[You can go here to view the search results.][hiv-search].
At the time of writing there were 302 results, but there are only a few that do not have any drug-like ligands: [`1TW7`](https://www.rcsb.org/structure/1TW7), [`2PC0`](https://www.rcsb.org/structure/2PC0), and [`2G69`](https://www.rcsb.org/structure/2G69).
Both [`1TW7`](https://www.rcsb.org/structure/1TW7) and [`2G69`](https://www.rcsb.org/structure/2G69) have one or more mutations studying drug resistance mechanism, so we will ignore these for now.
This leaves [`2PC0`](https://www.rcsb.org/structure/2PC0) as our protein.

<div id="hiv-protease-view" class="mol-container"></div>
<script>
var uri = './files/structures/2PC0.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#hiv-protease-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModelsAsFrames(data, "pdb");
        viewer.setStyle({}, {cartoon: {color: 'spectrum'}});
        viewer.center({chain: 'A'})
        viewer.setClickable({}, true, function(atom,viewer,event,container) {
            console.log(viewer.getView());
        });
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>

<!-- LINKS -->

[hiv-search]: https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22full_text%22%2C%22parameters%22%3A%7B%22value%22%3A%22HIV%20Protease%22%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%5D%2C%22logical_operator%22%3A%22and%22%2C%22label%22%3A%22full_text%22%7D%2C%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22entity_poly.rcsb_entity_polymer_type%22%2C%22value%22%3A%22Protein%22%2C%22operator%22%3A%22exact_match%22%7D%7D%5D%2C%22logical_operator%22%3A%22or%22%2C%22label%22%3A%22entity_poly.rcsb_entity_polymer_type%22%7D%2C%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_entry_info.resolution_combined%22%2C%22value%22%3A%7B%22from%22%3A0.5%2C%22to%22%3A2%2C%22include_lower%22%3Atrue%2C%22include_upper%22%3Atrue%7D%2C%22operator%22%3A%22range%22%7D%7D%5D%2C%22logical_operator%22%3A%22or%22%2C%22label%22%3A%22rcsb_entry_info.resolution_combined%22%7D%2C%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_polymer_entity.rcsb_ec_lineage.name%22%2C%22value%22%3A%22Hydrolases%22%2C%22operator%22%3A%22exact_match%22%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%5D%2C%22logical_operator%22%3A%22or%22%2C%22label%22%3A%22rcsb_polymer_entity.rcsb_ec_lineage.name%22%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%2C%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_entity_source_organism.ncbi_scientific_name%22%2C%22value%22%3A%22Human%20immunodeficiency%20virus%201%22%2C%22operator%22%3A%22exact_match%22%7D%7D%5D%2C%22logical_operator%22%3A%22or%22%2C%22label%22%3A%22rcsb_entity_source_organism.ncbi_scientific_name%22%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%5D%2C%22logical_operator%22%3A%22and%22%2C%22label%22%3A%22text%22%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%2C%22return_type%22%3A%22entry%22%2C%22request_options%22%3A%7B%22paginate%22%3A%7B%22start%22%3A0%2C%22rows%22%3A25%7D%2C%22results_content_type%22%3A%5B%22experimental%22%5D%2C%22sort%22%3A%5B%7B%22sort_by%22%3A%22score%22%2C%22direction%22%3A%22desc%22%7D%5D%2C%22scoring_strategy%22%3A%22combined%22%7D%2C%22request_info%22%3A%7B%22query_id%22%3A%229bf9fe45638e11d563f8a0807116b29a%22%7D%7D
