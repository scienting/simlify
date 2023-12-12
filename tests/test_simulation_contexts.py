from simlify.simulation.packages.amber.contexts import AmberContextValidator


class TestValidateContexts:
    def test_amber_protein_standard_context(self, amber_protein_standard_context):
        """Test if the standard amber protein simulation context is valid."""
        AmberContextValidator.validate(amber_protein_standard_context)
