# Marimo Notebook Example

## What is Marimo?

[Marimo](https://marimo.io/) is a reactive Python notebook that runs as an interactive web application. Unlike traditional Jupyter notebooks, Marimo notebooks are:

- **Reactive**: Cells automatically re-run when their dependencies change
- **Reproducible**: Notebooks are stored as pure Python files
- **Interactive**: Built-in UI elements for creating interactive apps

## Example Marimo Notebook

Below is an example of a Marimo notebook for computational chemistry. This would be saved as a `.py` file and can be run with `marimo edit filename.py`.

### Sample Code: Molecular Similarity Analysis

```python
import marimo

__generated_with = "0.1.0"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    return mo,


@app.cell
def __(mo):
    mo.md(
        """
        # Molecular Similarity Analysis
        
        This interactive notebook demonstrates how to calculate 
        molecular similarity using Tanimoto coefficients.
        """
    )
    return


@app.cell
def __():
    # Example SMILES strings
    molecules = {
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "Paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
    }
    return molecules,


@app.cell
def __(mo, molecules):
    # Create an interactive selector for molecules
    molecule_selector = mo.ui.dropdown(
        options=list(molecules.keys()),
        value="Aspirin",
        label="Select a molecule:"
    )
    
    mo.md(f"### Selected Molecule: {molecule_selector.value}")
    return molecule_selector,


@app.cell
def __(molecules, molecule_selector):
    selected_smiles = molecules[molecule_selector.value]
    
    # Display SMILES string
    print(f"SMILES: {selected_smiles}")
    return selected_smiles,


@app.cell
def __(mo):
    mo.md(
        """
        ## Calculating Molecular Properties
        
        In a full implementation, you would use RDKit to:
        1. Generate molecular fingerprints
        2. Calculate Tanimoto similarity
        3. Visualize molecular structures
        4. Compute physicochemical properties
        """
    )
    return


@app.cell
def __():
    # Conceptual example of similarity calculation
    def calculate_similarity(mol1_smiles, mol2_smiles):
        # In real implementation:
        # from rdkit import Chem
        # from rdkit.Chem import AllChem
        # from rdkit import DataStructs
        #
        # mol1 = Chem.MolFromSmiles(mol1_smiles)
        # mol2 = Chem.MolFromSmiles(mol2_smiles)
        #
        # fp1 = AllChem.GetMorganFingerprint(mol1, 2)
        # fp2 = AllChem.GetMorganFingerprint(mol2, 2)
        #
        # similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        # return similarity
        
        return 0.75  # Placeholder
    return calculate_similarity,


@app.cell
def __(mo):
    mo.md(
        """
        ## Benefits of Marimo for Computational Chemistry
        
        - **Reactive execution**: Changes propagate automatically
        - **No hidden state**: Reproducible results every time  
        - **Interactive widgets**: Easy to create parameter sweeps
        - **Pure Python files**: Easy to version control and share
        """
    )
    return


if __name__ == "__main__":
    app.run()
```

## Running Marimo Notebooks

To use Marimo notebooks:

1. **Install Marimo**:
   ```bash
   pip install marimo
   ```

2. **Create a new notebook**:
   ```bash
   marimo edit my_notebook.py
   ```

3. **Run an existing notebook**:
   ```bash
   marimo run my_notebook.py
   ```

## Converting to Static Content

For this blog, Marimo notebooks can be:

1. **Exported as HTML**: Use `marimo export html notebook.py > notebook.html`
2. **Included as code**: Show the Python source code (as above)
3. **Linked externally**: Host interactive versions separately

## Example Use Cases in Computational Chemistry

- **Parameter exploration**: Interactive widgets for exploring molecular properties
- **Data visualization**: Real-time plots of molecular descriptors
- **Structure analysis**: Interactive molecular similarity searches
- **QSAR modeling**: Interactive model building and validation

## Resources

- [Marimo Documentation](https://docs.marimo.io/)
- [Marimo GitHub](https://github.com/marimo-team/marimo)
- [Examples Gallery](https://marimo.io/gallery)

---

*Note: To see a live Marimo notebook, you need to run it locally or deploy it as a web application.*
