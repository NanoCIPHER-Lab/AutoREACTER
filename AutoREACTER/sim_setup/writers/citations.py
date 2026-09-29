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

    def write_citations(self):
        """Writes the collected citations to a .bib file."""
        filepath = Path(self.save_path) / "citations.bib"

        with open(filepath, "w") as f:
            # Join the citations and strip extra whitespace to format cleanly
            combined_citations = "\n\n".join(
                [c.strip() for c in self.citations]
            )
            f.write(combined_citations)
            f.write("\n")

        print(
            f"Successfully wrote {len(self.citations)} citations to {filepath}"
        )

    def get_citations_string(self) -> str:
        """Returns the compiled citations as a single string."""
        return "\n\n".join([c.strip() for c in self.citations])

    def cite_autoreacter(self):
        cite_text = """
@Software{AutoREACTER,
  author    = {Mahanthe, Janitha and Zulueta Nieto, Alessandra and Gissinger, Jacob},
  title     = {AutoREACTER},
  version   = {1.0.1},
  publisher = {Zenodo},
  year      = {2026},
  doi       = {10.5281/zenodo.23046173},
  url       = {https://doi.org/10.5281/zenodo.23046173},
}
"""
        return cite_text

    def cite_lunar(self):
        cite_text = """
@Article{Kemppainen24,
  author  = {Kemppainen, Josh and Gissinger, Jacob R. and Gowtham, S. and Odegard, Gregory M.},
  title   = {LUNAR: Automated Input Generation and Analysis for Reactive LAMMPS Simulations},
  journal = {Journal of Chemical Information and Modeling},
  year    = {2024},
  volume  = {64},
  number  = {13},
  pages   = {5108--5126},
  doi     = {10.1021/acs.jcim.4c00730},
}
"""
        return cite_text

    def cite_reacter(self):
        cite_text = """
@Article{Gissinger24cpc,
  author  = {Gissinger, J. R. and Jensen, B. D. and Wise, K. E.},
  title   = {Molecular Modeling of Reactive Systems with REACTER},
  journal = {Computer Physics Communications},
  year    = {2024},
  volume  = {304},
  pages   = {109287},
  doi     = {10.1016/j.cpc.2024.109287},
}
"""
        return cite_text