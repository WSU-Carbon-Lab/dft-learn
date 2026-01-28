"""Process test data from the output directory and save it as parquet files."""

from pathlib import Path

from dftlearn.io.stobe import process_directory


def main():
    """Process the test data and save it as parquet files."""
    test_data_dir = Path("tests/output")
    output_dir = Path("tests/reference_data")
    output_dir.mkdir(exist_ok=True)

    (
        basis_sets,
        energy_results,
        orbital_alpha,
        orbital_beta,
        xray_transitions,
        atomic_coordinates,
    ) = process_directory(
        directory=str(test_data_dir),
        width1=0.5,
        width2=1.0,
        ewid1=280,
        ewid2=300,
        verbose=False,
    )

    # Save outputs as parquet files
    basis_sets.to_parquet(output_dir / "basis_sets.parquet")
    energy_results.to_parquet(output_dir / "energy_results.parquet")
    orbital_alpha.to_parquet(output_dir / "orbital_alpha.parquet")
    orbital_beta.to_parquet(output_dir / "orbital_beta.parquet")
    xray_transitions.to_parquet(output_dir / "xray_transitions.parquet")
    atomic_coordinates.to_parquet(output_dir / "atomic_coordinates.parquet")


if __name__ == "__main__":
    main()
