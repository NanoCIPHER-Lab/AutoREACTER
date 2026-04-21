import os
from typing import TYPE_CHECKING

# Only import LunarFiles for type checking to avoid circular/runtime import overhead.
if TYPE_CHECKING:
    from AutoREACTER.reaction_preparation.lunar_client.lunar_api_wrapper import LunarFiles


class ForceFieldValidationError(Exception):
    """Raised when the force field file is missing, incomplete, or contains invalid coefficients."""
    pass


class FFValidator:
    """
    Validate a force field file used by LAMMPS-style workflows.

    The validator checks that:
    - the file exists,
    - required section labels are present,
    - each coefficient section contains data,
    - coefficient rows contain the expected number of numeric values,
    - no coefficient row is entirely zero.
    """

    def __init__(self, lunar_files: "LunarFiles"):
        """Store the file reference and validate the force field immediately."""
        self.lunar_files = lunar_files
        self.ff_file = self.lunar_files.force_field_data
        self.validate()

    def validate(self) -> bool:
        """
        Validate the force field file and raise an exception if any rule is violated.

        Returns:
            bool: True if validation succeeds.

        Raises:
            ForceFieldValidationError: If the file is missing or fails any validation check.
        """
        if not os.path.exists(self.ff_file):
            raise ForceFieldValidationError(f"File not found: {self.ff_file}")

        with open(self.ff_file, "r") as f:
            lines = f.readlines()

        # Define the sections we care about
        coeff_sections = ["Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs"]

        for section_name in coeff_sections:
            start_idx = self.find_section_start(lines, section_name)
            if start_idx is None:
                continue 

            data_start, data_end = self.find_data_block(lines, start_idx)
           
            self.check_coefficients(
                start_line=data_start,
                end_line=data_end,
                lines=lines,
                section_name=section_name,
            )
        return True

    def find_section_start(self, lines: list[str], section_name: str) -> int | None:
        """
        Find the first line containing a given section header.

        Args:
            lines: File contents split into individual lines.
            section_name: The section header to search for.

        Returns:
            The zero-based line index where the section starts, or None if not found.
        """
        for i, line in enumerate(lines):
            if section_name in line:
                return i
        return None

    def find_data_block(self, lines: list[str], section_start: int) -> tuple[int, int]:
        """
        Locate the coefficient data block below a section header.

        The method skips blank lines immediately after the header and then reads
        until the next blank line, which marks the end of the block.

        Example layout handled by this method:

            Bond Coeffs  # class2

            1 ...
            2 ...

        Args:
            lines: File contents split into individual lines.
            section_start: Line index of the section header.

        Returns:
            A tuple of (data_start, data_end), where data_end is exclusive.
        """
        i = section_start + 1

        # Skip blank lines directly after the section header.
        while i < len(lines) and not lines[i].strip():
            i += 1

        data_start = i

        # Read until the next blank line, which ends the data block.
        while i < len(lines) and lines[i].strip():
            i += 1

        data_end = i
        return data_start, data_end

    def check_coefficients(
        self,
        start_line: int,
        end_line: int,
        lines: list[str],
        section_name: str,
    ) -> bool:
        """
        Validate coefficient rows within a section.

        Each data row is expected to contain:
        - a type identifier,
        - the required number of coefficients,
        - numeric values only,
        - at least one coefficient that is non-zero.

        Args:
            start_line: First line of the data block.
            end_line: One past the last line of the data block.
            lines: File contents split into individual lines.
            section_name: Name of the section being validated.

        Returns:
            True if all rows in the section are valid.

        Raises:
            ForceFieldValidationError: If a row is malformed, non-numeric, or all zeros.
        """
        found_data = False
        expected_count = None  # This will be determined dynamically

        for line_no, line in enumerate(lines[start_line:end_line], start=start_line):
            # 1. Clean data (remove comments and whitespace)
            data = line.split("#")[0].strip()
            if not data:
                continue

            values = data.split()
            # values[0] is the ID (e.g., '1'), values[1:] are the coefficients
            current_coeffs = values[1:]

            # 2. Dynamic Count Logic
            if expected_count is None:
                expected_count = len(current_coeffs)
                # print(f"Section '{section_name}' dynamic count set to: {expected_count}")
            
            # 3. Validation against the dynamic count
            if len(current_coeffs) != expected_count:
                raise ForceFieldValidationError(
                    f"Inconsistent data in '{section_name}' at line {line_no+1}. "
                    f"Expected {expected_count} coefficients, found {len(current_coeffs)}."
                )

            # 4. Numeric and Zero-Check
            try:
                numeric_vals = [float(v) for v in current_coeffs]
                if all(v == 0.0 for v in numeric_vals):
                    raise ForceFieldValidationError(f"\nAll zeros in line {line_no+1} of {section_name} : \n\t{line.strip()}")
            except ValueError:
                raise ForceFieldValidationError(f"Non-numeric data in {section_name} at line {line_no+1}")

            found_data = True

        if not found_data:
            raise ForceFieldValidationError(f"No data found in section {section_name}")

        return True
