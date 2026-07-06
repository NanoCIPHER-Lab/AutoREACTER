MAX_LOOP = 5 # Maximum number of iterations for the reaction progression loop. Users should be able to adjust this value based on their specific needs and the complexity of the reactions being analyzed.
from typing import TYPE_CHECKING
from dataclasses import dataclass
if TYPE_CHECKING:
    from AutoREACTER.session import Session

from AutoREACTER.reaction_preparation.reaction_processor.detected_chemistry_filter import DetectedChemistryFilter, DetectedChemistries


@dataclass(slots=True)
class ReactionProgressionSession:
    """
    Holds the state of the reaction progression process, including detected chemistries and the current iteration count.
    """
    iteration: int = 0  # Current iteration count of the reaction progression loop.
    detected_chemistries: DetectedChemistries


class ReactionProgression:

    def __init__(self, session: "Session"):
        self.detected_chemistry_filter = DetectedChemistryFilter(session)
    def reaction_progression(self, session: "Session", max_loop: int = MAX_LOOP) -> None:
        pass