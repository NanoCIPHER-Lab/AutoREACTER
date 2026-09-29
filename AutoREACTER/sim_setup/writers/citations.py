from pathlib import Path
class CitationWriter:
    def __init__(self, save_path: str, LUNAR: bool = True):
        self.citations = []
        self.LUNAR = LUNAR
        self.save_path = save_path

        # Always include base citations
        self.citations.append(self.cite_autoreacter())
        self.citations.append(self.cite_reacter())
        
        # Include LUNAR citation if flag is True
        if self.LUNAR:
            self.citations.append(self.cite_lunar())

    def write_citations(self, ):
        """Writes the collected citations to a .bib file."""
        filepath: str = "citations.bib"
        filepath = Path(self.save_path) / filepath
        with open(filepath, 'w') as f:
            # Join the citations and strip extra whitespace to format cleanly
            combined_citations = "\n".join([c.strip() for c in self.citations])
            f.write(combined_citations)
            f.write("\n")
        print(f"Successfully wrote {len(self.citations)} citations to {filepath}")

    def get_citations_string(self) -> str:
        """Returns the compiled citations as a single string."""
        return "\n\n".join([c.strip() for c in self.citations])

    def cite_autoreacter(self):
        cite_text = """
@Misc{AutoREACTER,
  author       = {{NanoCIPHER Lab}},
  title        = {AutoREACTER},
  howpublished = {\\url{https://github.com/NanoCIPHER-Lab/AutoREACTER}},
  note         = {Citation pending -- DOI to be assigned via Zenodo on first release},
}
"""
        return cite_text

    def cite_lunar(self):
        cite_text = """
@Article{Kemppainen24,
  author  = {Kemppainen, Josh and Gissinger, Jacob R. and Gowtham, S. and Odegard, Gregory M.},
  title   = {LUNAR: Automated Input Generation and Analysis for Reactive LAMMPS Simulations},
  journal = {Journal of Chemical Information and Modeling},
  year    = 2024,
  volume  = 64,
  number  = 13,
  pages   = {5108--5126},
  doi     = {10.1021/acs.jcim.4c00730},
}
""" 
        return cite_text

    def cite_reacter(self):
        cite_text = """

@Article{Gissinger24cpc,
  author  = {J. R. Gissinger and B. D. Jensen and K. E. Wise},
  title   = {Molecular Modeling of Reactive Systems with REACTER},
  journal = {Computer Physics Communications},
  year    = 2024,
  volume  = 304,
  number  = 109287,
  doi     = {10.1016/j.cpc.2024.109287},
}
"""
        return cite_text