"""
XYZ Utilities for Molecular Structure Generation and Manipulation

This module provides tools for:
- Converting SMILES to XYZ coordinates (2D and 3D)
- Aligning molecules to specific axes
- Interactive relabeling of atoms
- Rotating molecular substructures
- Saving to XYZ file format
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# RDKit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.linear_model import LinearRegression

# Jupyter widget imports (optional - only needed for interactive tools)
try:
    import ipywidgets as widgets
    from IPython.display import clear_output, display

    WIDGETS_AVAILABLE = True
except ImportError:
    WIDGETS_AVAILABLE = False


# ==================== CONSTANTS ====================

COVALENT_RADII = {
    "C": 0.76,
    "H": 0.31,
    "S": 1.05,
    "O": 0.66,
    "N": 0.71,
    "P": 1.07,
    "F": 0.57,
    "Cl": 0.99,
    "Br": 1.14,
    "I": 1.33,
}

ATOMIC_COLORS = {
    "C": "black",
    "H": "lightgray",
    "O": "red",
    "N": "blue",
    "S": "yellow",
    "P": "orange",
    "F": "green",
    "Cl": "green",
    "Br": "brown",
    "I": "purple",
}


# ==================== UTILITY FUNCTIONS ====================


def get_atom_type(atom_label):
    """Extract atom type from label like 'C01' -> 'C'"""
    return "".join(c for c in atom_label if c.isalpha())


def load_xyz(filename):
    """
    Load an XYZ file and return a DataFrame.

    Parameters
    ----------
    filename : str
        Path to XYZ file

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ['atom', 'x', 'y', 'z']
    """
    df = pd.read_csv(filename, sep=r"\s+", header=None, skiprows=2)
    df.columns = ["atom", "x", "y", "z"]
    return df


def save_xyz_file(df, filename, comment="Generated XYZ structure"):
    """
    Save DataFrame to standard XYZ file format.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with columns ['atom', 'x', 'y', 'z']
    filename : str
        Output filename
    comment : str, optional
        Comment line in XYZ file
    """
    n_atoms = len(df)
    with open(filename, "w") as f:
        f.write(f"{n_atoms}\n")
        f.write(f"{comment}\n")
        for _, row in df.iterrows():
            x = np.float32(row["x"])
            y = np.float32(row["y"])
            z = np.float32(row["z"])
            f.write(f"{row['atom']:<4} {x:>12.6f} {y:>12.6f} {z:>12.6f}\n")
    print(f"✅ Saved {n_atoms} atoms to {filename}")


# ==================== SMILES TO XYZ CONVERSION ====================


def smiles_to_mol(smiles, add_hydrogens=True):
    """
    Convert SMILES string to RDKit molecule.

    Parameters
    ----------
    smiles : str
        SMILES string
    add_hydrogens : bool, optional
        Whether to add explicit hydrogens

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        RDKit molecule object
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    if add_hydrogens:
        mol = Chem.AddHs(mol)

    return mol


def rdkit_mol_to_xyz_2d(mol):
    """
    Convert RDKit molecule to XYZ coordinates in 2D (xy-plane).

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ['atom', 'x', 'y', 'z']
    """
    # Generate 2D coordinates
    AllChem.Compute2DCoords(mol)
    AllChem.MMFFOptimizeMolecule(mol, confId=0)
    conf = mol.GetConformer()

    # Extract atom information and coordinates
    atoms_data = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atom_symbol = atom.GetSymbol()

        atoms_data.append(
            {
                "atom": atom_symbol,
                "x": float(pos.x),
                "y": float(pos.y),
                "z": float(pos.z),
            }
        )

    # Create DataFrame
    xyz_df = pd.DataFrame(atoms_data)

    # Group atoms by type and add sequential numbering
    grouped_dfs = []
    for atom_type in sorted(xyz_df["atom"].unique()):
        atom_group = xyz_df[xyz_df["atom"] == atom_type].copy().reset_index(drop=True)
        n_atoms = len(atom_group)
        atom_group["atom"] = [f"{atom_type}{i + 1:02d}" for i in range(n_atoms)]
        grouped_dfs.append(atom_group)

    final_df = pd.concat(grouped_dfs, ignore_index=True)
    return final_df


def rdkit_mol_to_xyz_3d(mol):
    """
    Convert RDKit molecule to XYZ coordinates in 3D.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ['atom', 'x', 'y', 'z']
    """
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol, confId=0)
    conf = mol.GetConformer()

    # Extract atom information and coordinates
    atoms_data = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atom_symbol = atom.GetSymbol()

        atoms_data.append(
            {
                "atom": atom_symbol,
                "x": float(pos.x),
                "y": float(pos.y),
                "z": float(pos.z),
            }
        )

    # Create DataFrame
    xyz_df = pd.DataFrame(atoms_data)

    # Group atoms by type and add sequential numbering
    grouped_dfs = []
    for atom_type in sorted(xyz_df["atom"].unique()):
        atom_group = xyz_df[xyz_df["atom"] == atom_type].copy().reset_index(drop=True)
        n_atoms = len(atom_group)
        atom_group["atom"] = [f"{atom_type}{i + 1:02d}" for i in range(n_atoms)]
        grouped_dfs.append(atom_group)

    final_df = pd.concat(grouped_dfs, ignore_index=True)
    return final_df


def smiles_to_xyz(smiles, dimensions=2, add_hydrogens=True):
    """
    Convert SMILES string directly to XYZ DataFrame.

    Parameters
    ----------
    smiles : str
        SMILES string
    dimensions : int, optional
        2 for 2D coordinates, 3 for 3D coordinates
    add_hydrogens : bool, optional
        Whether to add explicit hydrogens

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ['atom', 'x', 'y', 'z']
    """
    mol = smiles_to_mol(smiles, add_hydrogens=add_hydrogens)

    if dimensions == 2:
        return rdkit_mol_to_xyz_2d(mol)
    elif dimensions == 3:
        return rdkit_mol_to_xyz_3d(mol)
    else:
        raise ValueError("dimensions must be 2 or 3")


# ==================== MOLECULE ALIGNMENT ====================


def center_molecule(df):
    """
    Center molecule at origin by subtracting centroid.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with columns ['atom', 'x', 'y', 'z']

    Returns
    -------
    pandas.DataFrame
        Centered DataFrame
    """
    df_centered = df.copy()
    coords = df_centered[["x", "y", "z"]].to_numpy()
    centroid = coords.mean(axis=0)
    df_centered[["x", "y", "z"]] = coords - centroid
    return df_centered


def rotate_about_axis(points, axis_point, angle):
    """
    Rotate points about an axis defined by axis_point and origin.

    Parameters
    ----------
    points : numpy.ndarray
        Array of 3D points to rotate
    axis_point : numpy.ndarray
        Point defining the rotation axis (from origin)
    angle : float
        Rotation angle in radians

    Returns
    -------
    numpy.ndarray
        Rotated points
    """
    # Translate points so axis_point is at the origin
    translated_points = points - axis_point

    # Normalize the rotation axis
    axis = axis_point / np.linalg.norm(axis_point)

    # Create the rotation matrix using Rodrigues' rotation formula
    K = np.array(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ]
    )
    I = np.identity(3)
    R = I + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

    # Apply the rotation matrix
    rotated_translated_points = np.dot(translated_points, R.T)

    # Translate points back
    rotated_points = rotated_translated_points + axis_point

    return rotated_points


def align_molecule(df, method="bond", **kwargs):
    """
    Align molecule by rotating to place specified features along the x-axis.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with columns ['atom', 'x', 'y', 'z']
    method : str
        'bond' - align a specific bond to x-axis
        'line' - fit line through specified atoms and align to x-axis
    **kwargs : dict
        For method='bond': atom1, atom2 (atom labels)
        For method='line': atoms (list of atom labels) or atom_type (str)

    Returns
    -------
    pandas.DataFrame
        DataFrame with aligned coordinates
    float
        Rotation angle applied (in degrees)
    """
    df_aligned = df.copy()

    if method == "bond":
        atom1 = kwargs.get("atom1")
        atom2 = kwargs.get("atom2")

        if not atom1 or not atom2:
            raise ValueError("For method='bond', specify both atom1 and atom2")

        # Get coordinates of the two atoms
        atom1_coords = df.loc[df["atom"] == atom1, ["x", "y"]].values
        atom2_coords = df.loc[df["atom"] == atom2, ["x", "y"]].values

        if len(atom1_coords) == 0 or len(atom2_coords) == 0:
            raise ValueError(f"Atoms {atom1} or {atom2} not found in molecule")

        # Compute the vector from atom1 to atom2
        vector = atom2_coords[0] - atom1_coords[0]

        # Compute rotation angle to align with x-axis
        angle = np.arctan2(vector[1], vector[0])

        # Translate to place atom1 at origin
        translated_coords = df_aligned[["x", "y"]].values - atom1_coords[0]

        print(f"Aligning {atom1}-{atom2} bond to x-axis")

    elif method == "line":
        if "atoms" in kwargs:
            atom_labels = kwargs["atoms"]
            mask = df["atom"].isin(atom_labels)
        elif "atom_type" in kwargs:
            atom_type = kwargs["atom_type"]
            if isinstance(atom_type, list):
                mask = df["atom"].str[0].isin(atom_type)
            else:
                mask = df["atom"].str.startswith(atom_type)
        else:
            raise ValueError("For method='line', specify either 'atoms' or 'atom_type'")

        if not mask.any():
            raise ValueError("No matching atoms found for line fitting")

        selected_coords = df.loc[mask, ["x", "y"]].values

        if len(selected_coords) < 2:
            raise ValueError("Need at least 2 atoms to fit a line")

        # Fit line
        model = LinearRegression()
        model.fit(selected_coords[:, 0].reshape(-1, 1), selected_coords[:, 1])
        slope = model.coef_[0]
        intercept = model.intercept_

        # Calculate angle to rotate line to x-axis
        angle = -np.arctan(slope)

        # Translate coordinates
        translated_coords = df_aligned[["x", "y"]].values.copy()
        translated_coords[:, 1] -= intercept

        print(f"Fitting line through {mask.sum()} atoms, slope: {slope:.4f}")

    else:
        raise ValueError("Method must be 'bond' or 'line'")

    # Create rotation matrix
    rotation_matrix = np.array(
        [[np.cos(-angle), -np.sin(-angle)], [np.sin(-angle), np.cos(-angle)]]
    )

    # Apply rotation
    rotated_coords = np.dot(translated_coords, rotation_matrix.T)

    # Update DataFrame
    df_aligned[["x", "y"]] = rotated_coords

    rotation_degrees = np.degrees(angle)
    print(f"Rotation angle: {rotation_degrees:.2f}°")

    return df_aligned, rotation_degrees


# ==================== VISUALIZATION ====================


def plot_molecule_2d(
    df,
    title="Molecule Structure",
    show_labels=True,
    hide_hydrogens=True,
    figsize=(10, 6),
):
    """
    Plot 2D molecule structure.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with columns ['atom', 'x', 'y', 'z']
    title : str, optional
        Plot title
    show_labels : bool, optional
        Whether to show atom labels
    hide_hydrogens : bool, optional
        Whether to hide hydrogen labels
    figsize : tuple, optional
        Figure size

    Returns
    -------
    matplotlib.figure.Figure, matplotlib.axes.Axes
        Figure and axes objects
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Get colors for atoms
    colors = [ATOMIC_COLORS.get(atom[0], "gray") for atom in df["atom"]]

    # Plot atoms
    ax.scatter(df["x"], df["y"], c=colors, s=300, edgecolors="black", alpha=0.8)

    # Add atom labels
    if show_labels:
        for i, row in df.iterrows():
            if hide_hydrogens and row["atom"].startswith("H"):
                continue
            label_color = "white" if row["atom"][0] in ["C", "N"] else "black"
            ax.annotate(
                row["atom"],
                (row["x"], row["y"]),
                ha="center",
                va="center",
                fontweight="bold",
                color=label_color,
                fontsize=8,
            )

    ax.set_xlabel("X coordinate (Å)")
    ax.set_ylabel("Y coordinate (Å)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")
    plt.tight_layout()

    return fig, ax


def calculate_bonds(df, tolerance=0.3):
    """
    Calculate bonds between atoms based on covalent radii.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with columns ['atom', 'x', 'y', 'z']
    tolerance : float, optional
        Tolerance for bond detection (Angstroms)

    Returns
    -------
    list of tuple
        List of (index1, index2) tuples representing bonds
    """
    bonds = []
    positions = df[["x", "y", "z"]].values

    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            idx1, idx2 = df.index[i], df.index[j]
            atom1_type = get_atom_type(df.loc[idx1, "atom"])
            atom2_type = get_atom_type(df.loc[idx2, "atom"])

            if atom1_type not in COVALENT_RADII or atom2_type not in COVALENT_RADII:
                continue

            expected_length = COVALENT_RADII[atom1_type] + COVALENT_RADII[atom2_type]
            actual_distance = np.linalg.norm(positions[i] - positions[j])

            if actual_distance <= expected_length + tolerance:
                bonds.append((idx1, idx2))

    return bonds


# ==================== INTERACTIVE RELABELER ====================


class MoleculeRelabeler:
    """
    Interactive widget for relabeling atoms in a molecule.
    Requires ipywidgets and matplotlib widget backend.
    """

    def __init__(self, df):
        if not WIDGETS_AVAILABLE:
            raise ImportError(
                "ipywidgets is required for MoleculeRelabeler. Install with: pip install ipywidgets"
            )

        self.original_df = df.copy().reset_index(drop=True)
        self._generate_initial_labels()

        self.df = self.original_df.copy()
        self.bonds = calculate_bonds(self.df)
        self.selected_atoms = []

        # Widget setup
        self.reset_button = widgets.Button(
            description="Reset All", button_style="warning"
        )
        self.swap_info = widgets.HTML(value="<b>Click atoms to select for swapping</b>")
        self.selected_atoms_display = widgets.HTML(value="Selected atoms: None")
        self.swap_button = widgets.Button(
            description="Swap Selected", button_style="info", disabled=True
        )
        self.clear_selection_button = widgets.Button(
            description="Clear Selection", button_style="primary"
        )
        self.status_output = widgets.Output()

        # Plotting setup
        plt.ioff()
        self.fig, self.ax = plt.subplots(figsize=(10, 7))
        plt.ion()

        self.ax.set_xlabel("X-coordinate (Å)")
        self.ax.set_ylabel("Y-coordinate (Å)")
        self.ax.set_title("Atom Swapping Tool")
        self.ax.grid(True, linestyle="--", alpha=0.3)
        self.ax.set_aspect("equal", adjustable="box")

        self.atom_scatter = self.ax.scatter([], [], edgecolors="k", s=300, zorder=3)
        self.highlight_scatter = self.ax.scatter(
            [], [], c="none", edgecolors="cyan", s=600, linewidths=3, zorder=4
        )
        self.bond_lines = []
        self.atom_annotations = []

        self.fig.tight_layout()

        # Event handlers
        self.reset_button.on_click(self.reset_labels)
        self.swap_button.on_click(self.swap_atoms)
        self.clear_selection_button.on_click(self.clear_selection)
        self.fig.canvas.mpl_connect("button_press_event", self.on_plot_click)

        # Layout
        controls = widgets.VBox(
            [
                widgets.HBox([self.reset_button, self.clear_selection_button]),
                widgets.HTML("<hr>"),
                self.swap_info,
                self.selected_atoms_display,
                self.swap_button,
            ]
        )
        controls.layout = widgets.Layout(width="30%", padding="10px")

        main_ui = widgets.HBox([controls, self.fig.canvas])
        final_ui = widgets.VBox([main_ui, self.status_output])

        display(final_ui)
        self.update_plot()

    def _generate_initial_labels(self):
        """Generate initial sequential labels for atoms."""
        atom_counts = {}
        new_labels = []
        for atom_symbol in self.original_df["atom"]:
            atom_type = get_atom_type(atom_symbol)
            atom_counts.setdefault(atom_type, 0)
            atom_counts[atom_type] += 1
            new_labels.append(f"{atom_type}{atom_counts[atom_type]:02d}")
        self.original_df["atom"] = new_labels

    def reset_labels(self, button):
        """Reset to original atom positions."""
        self.df = self.original_df.copy()
        self.bonds = calculate_bonds(self.df)
        with self.status_output:
            clear_output(wait=True)
            print("Molecule reset to original order.")
        self.clear_selection()

    def on_plot_click(self, event):
        """Handle plot click events."""
        if event.button != 1 or not event.inaxes:
            return

        plot_df = self.df[[get_atom_type(atom) != "H" for atom in self.df["atom"]]]
        if plot_df.empty:
            return

        click_x, click_y = event.xdata, event.ydata

        distances = []
        for idx, row in plot_df.iterrows():
            dist = np.sqrt((row["x"] - click_x) ** 2 + (row["y"] - click_y) ** 2)
            distances.append((dist, idx))

        if not distances:
            return
        distances.sort()
        closest_dist, closest_idx = distances[0]

        if closest_dist < 0.5:
            if closest_idx in self.selected_atoms:
                self.selected_atoms.remove(closest_idx)
            else:
                self.selected_atoms.append(closest_idx)
                if len(self.selected_atoms) > 2:
                    self.selected_atoms.pop(0)
            self.update_selection_display()
            self.update_plot()

    def update_selection_display(self):
        """Update the selection display widget."""
        if not self.selected_atoms:
            self.selected_atoms_display.value = "Selected atoms: None"
            self.swap_button.disabled = True
        elif len(self.selected_atoms) == 1:
            atom_label = self.df.loc[self.selected_atoms[0], "atom"]
            self.selected_atoms_display.value = (
                f"Selected: <b>{atom_label}</b> (select one more)"
            )
            self.swap_button.disabled = True
        else:
            atom1_label = self.df.loc[self.selected_atoms[0], "atom"]
            atom2_label = self.df.loc[self.selected_atoms[1], "atom"]
            self.selected_atoms_display.value = (
                f"Selected: <b>{atom1_label}</b> ↔ <b>{atom2_label}</b>"
            )
            self.swap_button.disabled = False

    def clear_selection(self, button=None):
        """Clear atom selection."""
        self.selected_atoms = []
        self.update_selection_display()
        self.update_plot()

    def swap_atoms(self, button):
        """Swap positions of selected atoms."""
        if len(self.selected_atoms) != 2:
            return
        idx1, idx2 = self.selected_atoms

        old_label1 = self.df.loc[idx1, "atom"]
        old_label2 = self.df.loc[idx2, "atom"]

        # Swap coordinates
        coords1 = self.df.loc[idx1, ["x", "y", "z"]].copy()
        coords2 = self.df.loc[idx2, ["x", "y", "z"]].copy()

        self.df.loc[idx1, ["x", "y", "z"]] = coords2
        self.df.loc[idx2, ["x", "y", "z"]] = coords1

        self.bonds = calculate_bonds(self.df)
        self.update_plot()

        with self.status_output:
            clear_output(wait=True)
            print(f"✅ Swapped positions of atoms: {old_label1} ↔ {old_label2}")
            print("    Labels remain the same, only coordinates were swapped")

        self.clear_selection()

    def update_plot(self, button=None):
        """Update the plot with current molecule state."""
        plot_df = self.df[[get_atom_type(atom) != "H" for atom in self.df["atom"]]]

        # Update atom scatter
        atom_positions = (
            plot_df[["x", "y"]].values if not plot_df.empty else np.empty((0, 2))
        )
        atom_types = [get_atom_type(atom) for atom in plot_df["atom"]]
        colors = [ATOMIC_COLORS.get(t, "pink") for t in atom_types]
        self.atom_scatter.set_offsets(atom_positions)
        self.atom_scatter.set_facecolors(colors)

        # Update bonds and annotations
        for line in self.bond_lines:
            line.remove()
        self.bond_lines.clear()

        for ann in self.atom_annotations:
            ann.remove()
        self.atom_annotations.clear()

        visible_indices = set(plot_df.index)
        for idx1, idx2 in self.bonds:
            if idx1 in visible_indices and idx2 in visible_indices:
                x_coords = [self.df.loc[idx1, "x"], self.df.loc[idx2, "x"]]
                y_coords = [self.df.loc[idx1, "y"], self.df.loc[idx2, "y"]]
                (line,) = self.ax.plot(
                    x_coords, y_coords, "grey", alpha=0.7, linewidth=4, zorder=1
                )
                self.bond_lines.append(line)

        for idx, row in plot_df.iterrows():
            label_color = (
                "white" if get_atom_type(row["atom"]) in ["C", "N"] else "black"
            )
            ann = self.ax.annotate(
                row["atom"],
                (row["x"], row["y"]),
                ha="center",
                va="center",
                fontsize=9,
                color=label_color,
                fontweight="bold",
                zorder=5,
            )
            self.atom_annotations.append(ann)

        # Update highlight
        if self.selected_atoms:
            selected_data = self.df.loc[self.selected_atoms, ["x", "y"]].values
            self.highlight_scatter.set_offsets(selected_data)
        else:
            self.highlight_scatter.set_offsets(np.empty((0, 2)))

        self.ax.relim()
        self.ax.autoscale_view()
        self.fig.canvas.draw_idle()

    def get_updated_dataframe(self):
        """Return the updated dataframe with current atom positions."""
        return self.df.copy()
