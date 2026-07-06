import unittest
#this should import the main class from the code
from AutoREACTER.input_parser import (
    InputParser,
    NumericFieldError,
    InputSchemaError,
    CompatibilityError,
    DuplicateMonomerError
)

class TestInputParser(unittest.TestCase) :
    
    def setUp(self):
        """Runs automatically before every test to give us a fresh parser."""
        self.parser = InputParser()
    #===============
    # Tests for:  _validate_temperature
    #===============

    def test_validate_temperature_valid(self):
        """Good Ending: A valid positive temp should return as a float."""
        result = self.parser._validate_temperature(300)
        self.assertEqual(result, 300.0)
    
    def test_validate_temperature_negative(self):
        """Bad Ending: A negative temp has to raise a NumericFieldError."""
        with self.assertRaises(NumericFieldError):
            self.parser._validate_temperature(-150)

    def test_validate_temperature_zero(self):
        """Bad Ending: Absolute zero has to raise a NumericFieldError."""
        with self.assertRaises(NumericFieldError):
            self.parser._validate_temperature(0)

    def test_validate_temperature_wrong_type(self):
        """Bad Ending: If a string or boolean is given instead of a number it has to fail."""
        with self.assertRaises(NumericFieldError):
            self.parser._validate_temperature("room_temperature")

    #===============
    # Tests for:  _validate_density
    #===============   
    def test_validate_density_valid(self):
        """Good Ending: Valid positive density returns as float."""
        result = self.parser._validate_density(0.85)
        self.assertEqual(result, 0.85)
    
    def test_validate_density_negative_or_zero(self):
        """Bad Ending: Rejects zero or negative density values."""
        with self.assertRaises(NumericFieldError):
            self.parser._validate_density(0)
        with self.assertRaises(NumericFieldError):
            self.parser._validate_density(-0.5)
    
    #================
    # Tests for: _validate_force_field
    #================

    def test_validate_force_field_default(self):
        """Good Ending: If None then defaults to 'PCFF'."""
        result = self.parser._validate_force_field(None)
        self.assertEqual(result, "PCFF")
    
    def test_validate_force_field_canonical(self):
        """Good path: Normalizes force-field capitalization and aliases."""
        result = self.parser._validate_force_field("pcff-iff")
        self.assertEqual(result, "PCFF-IFF")

    def test_validate_force_field_unsupported (self):
        """Bad Ending: Passing a completely random name should raise an InputSchemaError."""
        with self.assertRaises(InputSchemaError):
            self.parser._validate_force_field("NotAForceField")

    def test_validate_force_field_incompatible(self):
        """Bad Ending: 'OPLSAA' is recognized but incompatible, so it should raise a CompatibilityError."""
        with self.assertRaises(CompatibilityError):
            self.parser._validate_force_field("oplsaa")
    
    # ==========================================
    # Tests for: validate_no_duplicate_smiles
    # ==========================================

    def test_validate_no_duplicate_smiles_good_path(self):
        """Good Ending: Adding a unique SMILES appends it to the tracker list."""
        tracker_list = ["CCO", "C1=CC(=CC(=C1)N)N"]
        new_smiles = "ClC(=O)c1cc(cc(c1)C(Cl)=O)C(Cl)=O"
        
        result = self.parser.validate_no_duplicate_smiles(new_smiles, tracker_list)
        
        # Verify the list grew to 3 items and includes this new molecule
        self.assertEqual(len(result), 3)
        self.assertIn(new_smiles, result)

    def test_validate_no_duplicate_smiles_raises_error(self):
        """Bad Ending: Adding an existing SMILES triggers DuplicateMonomerError."""
        tracker_list = ["CCO", "C1=CC(=CC(=C1)N)N"]
        duplicate_smiles = "CCO"  # Already present!
        
        with self.assertRaises(DuplicateMonomerError):
            self.parser.validate_no_duplicate_smiles(duplicate_smiles, tracker_list)

if __name__ == "__main__":
    unittest.main() 
