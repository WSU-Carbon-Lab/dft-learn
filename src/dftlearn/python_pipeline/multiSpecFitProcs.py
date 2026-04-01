import numpy as np
import pandas as pd
from scipy.optimize import minimize, Bounds
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
import warnings
warnings.filterwarnings('ignore')

class MultiSpectrumGaussianFitter:
    """
    Simultaneous fitting of multiple spectra using Gaussian peaks with staged optimization.
    
    Stages:
    1. Optimize only amplitudes (OS)
    2. Optimize amplitudes and positions (E)
    3. Optimize amplitudes, positions, and widths
    """
    
    def __init__(self, parameters_file: str, experimental_file: str):
        """
        Initialize the fitter with Gaussian parameters and experimental data.
        
        Args:
            parameters_file: CSV file with Gaussian parameters
            experimental_file: CSV file with experimental spectra
        """
        self.load_parameters(parameters_file)
        self.load_experimental_data(experimental_file)
        self.setup_initial_parameters()
        
    def load_parameters(self, filename: str):
        """Load Gaussian parameters from CSV file."""
        self.params_df = pd.read_csv(filename)
        self.n_gaussians = len(self.params_df)
        
        # Extract relevant parameters
        self.initial_positions = self.params_df['E'].values
        self.initial_widths = self.params_df['width'].values
        self.initial_amplitudes = self.params_df['OS'].values
        
        print(f"Loaded {self.n_gaussians} Gaussian peaks")
        print(f"Energy range: {self.initial_positions.min():.2f} - {self.initial_positions.max():.2f} eV")
        
    def load_experimental_data(self, filename: str):
        """Load experimental spectra data."""
        self.exp_df = pd.read_csv(filename)
        
        # Identify energy and intensity columns
        self.energy_columns = [col for col in self.exp_df.columns if col.startswith('E_')]
        self.intensity_columns = [col for col in self.exp_df.columns if col.startswith('CuPc_CuI_')]
        
        self.angles = [col.split('_')[1] for col in self.energy_columns]
        self.n_spectra = len(self.angles)
        
        # Store experimental data for each spectrum
        self.exp_energies = {}
        self.exp_intensities = {}
        
        for i, angle in enumerate(self.angles):
            energy_col = self.energy_columns[i]
            intensity_col = self.intensity_columns[i]
            
            # Remove NaN values
            mask = ~(pd.isna(self.exp_df[energy_col]) | pd.isna(self.exp_df[intensity_col]))
            self.exp_energies[angle] = self.exp_df.loc[mask, energy_col].values
            self.exp_intensities[angle] = self.exp_df.loc[mask, intensity_col].values
            
        print(f"Loaded {self.n_spectra} experimental spectra at angles: {self.angles}")
        
    def setup_initial_parameters(self):
        """Set up initial parameter arrays for optimization."""
        # Parameter organization: [amplitudes, positions, widths] for each spectrum
        # Each spectrum has its own amplitude scaling, but shares positions and widths
        
        # For each spectrum, we have separate amplitudes
        self.n_params_per_spectrum = self.n_gaussians  # amplitudes only
        self.n_shared_params = 2 * self.n_gaussians  # positions + widths
        
        # Total parameters: amplitudes for each spectrum + shared positions + shared widths
        self.total_params = self.n_spectra * self.n_gaussians + self.n_shared_params
        
    def gaussian(self, x: np.ndarray, position: float, width: float, amplitude: float) -> np.ndarray:
        """
        Calculate Gaussian peak values.
        
        Args:
            x: Energy array
            position: Peak position (eV)
            width: Peak width (FWHM)
            amplitude: Peak amplitude (oscillator strength)
            
        Returns:
            Gaussian values
        """
        sigma = width / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to sigma
        return amplitude * np.exp(-0.5 * ((x - position) / sigma) ** 2)
        
    def calculate_spectrum(self, energies: np.ndarray, amplitudes: np.ndarray, 
                          positions: np.ndarray, widths: np.ndarray) -> np.ndarray:
        """
        Calculate complete spectrum from all Gaussian peaks.
        
        Args:
            energies: Energy array
            amplitudes: Amplitude array for all peaks
            positions: Position array for all peaks  
            widths: Width array for all peaks
            
        Returns:
            Calculated spectrum intensity
        """
        spectrum = np.zeros_like(energies)
        for i in range(self.n_gaussians):
            spectrum += self.gaussian(energies, positions[i], widths[i], amplitudes[i])
        return spectrum
        
    def unpack_parameters(self, params: np.ndarray, stage: str) -> Tuple[Dict, np.ndarray, np.ndarray]:
        """
        Unpack parameter array based on optimization stage.
        
        Args:
            params: Flattened parameter array
            stage: Optimization stage ('amplitudes', 'positions', 'widths')
            
        Returns:
            amplitudes_dict: Dictionary of amplitudes for each spectrum
            positions: Shared positions array
            widths: Shared widths array
        """
        amplitudes_dict = {}
        
        if stage == 'amplitudes':
            # Only amplitudes are being optimized
            for i, angle in enumerate(self.angles):
                start_idx = i * self.n_gaussians
                end_idx = (i + 1) * self.n_gaussians
                amplitudes_dict[angle] = params[start_idx:end_idx]
            positions = self.initial_positions.copy()
            widths = self.initial_widths.copy()
            
        elif stage == 'positions':
            # Amplitudes and positions are being optimized
            for i, angle in enumerate(self.angles):
                start_idx = i * self.n_gaussians
                end_idx = (i + 1) * self.n_gaussians
                amplitudes_dict[angle] = params[start_idx:end_idx]
            
            pos_start = self.n_spectra * self.n_gaussians
            pos_end = pos_start + self.n_gaussians
            positions = params[pos_start:pos_end]
            widths = self.initial_widths.copy()
            
        elif stage == 'widths':
            # All parameters are being optimized
            for i, angle in enumerate(self.angles):
                start_idx = i * self.n_gaussians
                end_idx = (i + 1) * self.n_gaussians
                amplitudes_dict[angle] = params[start_idx:end_idx]
                
            pos_start = self.n_spectra * self.n_gaussians
            pos_end = pos_start + self.n_gaussians
            positions = params[pos_start:pos_end]
            
            width_start = pos_end
            width_end = width_start + self.n_gaussians
            widths = params[width_start:width_end]
            
        return amplitudes_dict, positions, widths
        
    def objective_function(self, params: np.ndarray, stage: str) -> float:
        """
        Objective function for least squares fitting.
        
        Args:
            params: Parameter array
            stage: Optimization stage
            
        Returns:
            Sum of squared residuals
        """
        amplitudes_dict, positions, widths = self.unpack_parameters(params, stage)
        
        total_residual = 0.0
        
        for angle in self.angles:
            # Calculate theoretical spectrum
            calc_spectrum = self.calculate_spectrum(
                self.exp_energies[angle], 
                amplitudes_dict[angle], 
                positions, 
                widths
            )
            
            # Calculate residuals
            residuals = self.exp_intensities[angle] - calc_spectrum
            total_residual += np.sum(residuals ** 2)
            
        return total_residual
        
    def setup_bounds_and_constraints(self, stage: str) -> Bounds:
        """
        Set up parameter bounds based on optimization stage.
        
        Args:
            stage: Optimization stage
            
        Returns:
            Bounds object for scipy optimization
        """
        lower_bounds = []
        upper_bounds = []
        
        # Amplitude bounds (must be positive)
        for i in range(self.n_spectra):
            lower_bounds.extend([0.0] * self.n_gaussians)
            upper_bounds.extend([np.inf] * self.n_gaussians)
            
        if stage in ['positions', 'widths']:
            # Position bounds (initial ± 2*width)
            for i in range(self.n_gaussians):
                lower_bound = self.initial_positions[i] - 2 * self.initial_widths[i]
                upper_bound = self.initial_positions[i] + 2 * self.initial_widths[i]
                lower_bounds.append(lower_bound)
                upper_bounds.append(upper_bound)
                
        if stage == 'widths':
            # Width bounds (positive, not greater than 2*initial)
            for i in range(self.n_gaussians):
                lower_bounds.append(1e-6)  # Small positive value
                upper_bounds.append(2 * self.initial_widths[i])
                
        return Bounds(lower_bounds, upper_bounds)
        
    def get_initial_parameters(self, stage: str) -> np.ndarray:
        """
        Get initial parameter array for optimization stage.
        
        Args:
            stage: Optimization stage
            
        Returns:
            Initial parameter array
        """
        params = []
        
        # Add amplitudes for each spectrum (start with initial values)
        for i in range(self.n_spectra):
            params.extend(self.initial_amplitudes)
            
        if stage in ['positions', 'widths']:
            # Add positions
            params.extend(self.initial_positions)
            
        if stage == 'widths':
            # Add widths
            params.extend(self.initial_widths)
            
        return np.array(params)
        
    def fit_stage(self, stage: str, initial_params: np.ndarray = None) -> Dict:
        """
        Perform fitting for a specific stage.
        
        Args:
            stage: Optimization stage
            initial_params: Initial parameters (if None, use default)
            
        Returns:
            Dictionary with optimization results
        """
        print(f"\n=== Stage: {stage.upper()} ===")
        
        if initial_params is None:
            initial_params = self.get_initial_parameters(stage)
            
        bounds = self.setup_bounds_and_constraints(stage)
        
        print(f"Optimizing {len(initial_params)} parameters...")
        print(f"Initial objective: {self.objective_function(initial_params, stage):.6f}")
        
        # Perform optimization
        result = minimize(
            fun=self.objective_function,
            x0=initial_params,
            args=(stage,),
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': 1000, 'ftol': 1e-9}
        )
        
        print(f"Final objective: {result.fun:.6f}")
        print(f"Optimization success: {result.success}")
        print(f"Number of iterations: {result.nit}")
        
        return {
            'result': result,
            'stage': stage,
            'final_params': result.x,
            'objective_value': result.fun,
            'success': result.success
        }
        
    def fit_all_stages(self) -> Dict:
        """
        Perform complete staged fitting: amplitudes -> positions -> widths.
        
        Returns:
            Dictionary with all fitting results
        """
        print("Starting staged optimization...")
        
        results = {}
        
        # Stage 1: Optimize amplitudes only
        stage1_result = self.fit_stage('amplitudes')
        results['stage1'] = stage1_result
        
        # Stage 2: Optimize amplitudes and positions
        # Use stage 1 results as starting point
        stage2_initial = self.get_initial_parameters('positions')
        # Update amplitudes from stage 1
        n_amp_params = self.n_spectra * self.n_gaussians
        stage2_initial[:n_amp_params] = stage1_result['final_params'][:n_amp_params]
        
        stage2_result = self.fit_stage('positions', stage2_initial)
        results['stage2'] = stage2_result
        
        # Stage 3: Optimize all parameters
        # Use stage 2 results as starting point
        stage3_initial = self.get_initial_parameters('widths')
        stage3_initial[:len(stage2_result['final_params'])] = stage2_result['final_params']
        
        stage3_result = self.fit_stage('widths', stage3_initial)
        results['stage3'] = stage3_result
        
        # Store final parameters
        self.final_parameters = stage3_result['final_params']
        self.final_amplitudes, self.final_positions, self.final_widths = \
            self.unpack_parameters(self.final_parameters, 'widths')
            
        return results
        
    def calculate_fit_statistics(self) -> Dict:
        """
        Calculate fitting statistics for final parameters.
        
        Returns:
            Dictionary with fit statistics
        """
        if not hasattr(self, 'final_parameters'):
            raise ValueError("Must run fitting first!")
            
        stats = {}
        total_points = 0
        total_ss_res = 0
        total_ss_tot = 0
        
        for angle in self.angles:
            # Calculate fitted spectrum
            fitted = self.calculate_spectrum(
                self.exp_energies[angle],
                self.final_amplitudes[angle],
                self.final_positions,
                self.final_widths
            )
            
            # Calculate statistics
            residuals = self.exp_intensities[angle] - fitted
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((self.exp_intensities[angle] - np.mean(self.exp_intensities[angle])) ** 2)
            
            n_points = len(self.exp_intensities[angle])
            rmse = np.sqrt(ss_res / n_points)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            stats[f'angle_{angle}'] = {
                'n_points': n_points,
                'rmse': rmse,
                'r_squared': r_squared,
                'ss_res': ss_res,
                'ss_tot': ss_tot
            }
            
            total_points += n_points
            total_ss_res += ss_res
            total_ss_tot += ss_tot
            
        # Overall statistics
        stats['overall'] = {
            'total_points': total_points,
            'total_rmse': np.sqrt(total_ss_res / total_points),
            'total_r_squared': 1 - (total_ss_res / total_ss_tot) if total_ss_tot > 0 else 0,
            'n_parameters': len(self.final_parameters)
        }
        
        return stats
    

import numpy as np
import pandas as pd
from scipy.optimize import minimize, Bounds
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Optional
import warnings
warnings.filterwarnings('ignore')

class NEXAFSTensorFitter:
    """
    Simultaneous fitting of multiple NEXAFS spectra using molecular absorption tensors.
    
    This implementation follows the IGOR code approach where each transition is modeled
    as a molecular absorption tensor that is transformed based on sample tilt (alpha)
    and measurement geometry (theta, phi).
    """
    
    def __init__(self, parameters_file: str, experimental_file: str):
        """
        Initialize the fitter with tensor parameters and experimental data.
        
        Args:
            parameters_file: CSV file with Gaussian parameters including tensor elements
            experimental_file: CSV file with experimental spectra at different angles
        """
        self.load_parameters(parameters_file)
        self.load_experimental_data(experimental_file)
        self.setup_tensor_parameters()
        
    def load_parameters(self, filename: str):
        """Load initial parameters from CSV file."""
        self.params_df = pd.read_csv(filename)
        self.n_transitions = len(self.params_df)
        
        # Extract basic parameters
        self.initial_positions = self.params_df['E'].values
        self.initial_widths = self.params_df['width'].values  
        self.initial_amplitudes = self.params_df['OS'].values
        
        # Extract tensor elements if available, otherwise use isotropic assumption
        tensor_cols = ['Txx', 'Txy', 'Txz', 'Tyy', 'Tyz', 'Tzz']
        if all(col in self.params_df.columns for col in tensor_cols):
            self.initial_tensor_elements = self.params_df[tensor_cols].values
        else:
            # Default to isotropic tensors with in-plane/out-of-plane ratio
            print("Tensor elements not found, using isotropic assumption")
            self.initial_tensor_elements = self._create_default_tensors()
            
        print(f"Loaded {self.n_transitions} transitions")
        print(f"Energy range: {self.initial_positions.min():.2f} - {self.initial_positions.max():.2f} eV")
        
    def _create_default_tensors(self):
        """Create default isotropic tensors with adjustable anisotropy."""
        tensors = np.zeros((self.n_transitions, 6))  # [Txx, Txy, Txz, Tyy, Tyz, Tzz]
        
        for i in range(self.n_transitions):
            # Assume Txx = Tyy (in-plane), Tzz (out-of-plane), others zero
            in_plane = self.initial_amplitudes[i] * 0.7  # Default 70% in-plane
            out_plane = self.initial_amplitudes[i] * 0.3  # Default 30% out-of-plane
            
            tensors[i] = [in_plane, 0, 0, in_plane, 0, out_plane]
            
        return tensors
        
    def load_experimental_data(self, filename: str):
        """Load experimental spectra data."""
        self.exp_df = pd.read_csv(filename)
        
        # Identify energy and intensity columns
        self.energy_columns = [col for col in self.exp_df.columns if col.startswith('E_')]
        self.intensity_columns = [col for col in self.exp_df.columns if col.startswith('CuPc_CuI_')]
        
        # Extract angles from column names
        self.angles = []
        for col in self.energy_columns:
            angle_str = col.split('_')[1]
            self.angles.append(float(angle_str))
            
        self.n_spectra = len(self.angles)
        
        # Store experimental data
        self.exp_energies = {}
        self.exp_intensities = {}
        
        for i, angle in enumerate(self.angles):
            energy_col = self.energy_columns[i]
            intensity_col = self.intensity_columns[i]
            
            # Remove NaN values
            mask = ~(pd.isna(self.exp_df[energy_col]) | pd.isna(self.exp_df[intensity_col]))
            self.exp_energies[angle] = self.exp_df.loc[mask, energy_col].values
            self.exp_intensities[angle] = self.exp_df.loc[mask, intensity_col].values
            
        print(f"Loaded {self.n_spectra} spectra at angles: {self.angles}")
        
    def setup_tensor_parameters(self):
        """Set up parameter organization for tensor-based fitting."""
        # Global parameters: i0 (intensity scale), alpha (tilt), phi (azimuth)
        self.n_global_params = 3
        
        # Per-transition parameters: position, width, tensor elements (6)
        self.n_params_per_transition = 8  # E, width, Txx, Txy, Txz, Tyy, Tyz, Tzz
        
        # Total parameters
        self.total_params = self.n_global_params + self.n_transitions * self.n_params_per_transition
        
        print(f"Total parameters to fit: {self.total_params}")
        
    def create_molecular_tensor(self, tensor_elements: np.ndarray) -> np.ndarray:
        """
        Create 3x3 molecular absorption tensor from 6 unique elements.
        
        Args:
            tensor_elements: [Txx, Txy, Txz, Tyy, Tyz, Tzz]
            
        Returns:
            3x3 symmetric tensor matrix
        """
        tensor = np.zeros((3, 3))
        tensor[0, 0] = tensor_elements[0]  # Txx
        tensor[0, 1] = tensor[1, 0] = tensor_elements[1]  # Txy
        tensor[0, 2] = tensor[2, 0] = tensor_elements[2]  # Txz
        tensor[1, 1] = tensor_elements[3]  # Tyy
        tensor[1, 2] = tensor[2, 1] = tensor_elements[4]  # Tyz
        tensor[2, 2] = tensor_elements[5]  # Tzz
        return tensor
        
    def create_rotation_matrices(self):
        """Create rotation matrices for azimuthal averaging."""
        angles = [0, 90, 180, 270]
        matrices = []
        
        for angle_deg in angles:
            angle_rad = np.radians(angle_deg)
            cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
            
            rot_matrix = np.array([
                [cos_a, sin_a, 0],
                [-sin_a, cos_a, 0],
                [0, 0, 1]
            ])
            matrices.append(rot_matrix)
            
        return matrices
        
    def create_tilt_matrix(self, alpha: float) -> np.ndarray:
        """
        Create tilt matrix for sample rotation around x-axis.
        
        Args:
            alpha: Tilt angle in radians
            
        Returns:
            3x3 rotation matrix
        """
        cos_a, sin_a = np.cos(alpha), np.sin(alpha)
        return np.array([
            [1, 0, 0],
            [0, cos_a, sin_a],
            [0, -sin_a, cos_a]
        ])
        
    def transform_molecular_tensor(self, molecular_tensor: np.ndarray, 
                                 alpha: float) -> np.ndarray:
        """
        Transform molecular tensor to film tensor including tilt and azimuthal averaging.
        
        Args:
            molecular_tensor: 3x3 molecular absorption tensor
            alpha: Sample tilt angle in radians
            
        Returns:
            3x3 film-averaged tensor
        """
        # Create tilt matrix
        tilt_matrix = self.create_tilt_matrix(alpha)
        
        # Apply tilt transformation
        tilted_tensor = tilt_matrix @ molecular_tensor @ tilt_matrix.T
        
        # Apply azimuthal averaging (4-fold rotation)
        rotation_matrices = self.create_rotation_matrices()
        film_tensor = np.zeros((3, 3))
        
        for rot_matrix in rotation_matrices:
            rotated_tensor = rot_matrix @ tilted_tensor @ rot_matrix.T
            film_tensor += rotated_tensor
            
        return film_tensor / 4.0  # Average over 4 rotations
        
    def calculate_mass_absorption(self, film_tensor: np.ndarray, 
                                theta: float, phi: float) -> float:
        """
        Calculate mass absorption for given measurement geometry.
        
        Args:
            film_tensor: 3x3 film-averaged tensor
            theta: Polar angle in radians
            phi: Azimuthal angle in radians
            
        Returns:
            Mass absorption value
        """
        # Electric field vector
        e_field = np.array([
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta)
        ])
        
        # Calculate absorption: E^T · T · E
        absorption = e_field.T @ film_tensor @ e_field
        return max(0, absorption)  # Ensure non-negative
        
    def gaussian_peak(self, energies: np.ndarray, position: float, 
                     width: float, amplitude: float) -> np.ndarray:
        """
        Calculate Gaussian peak profile.
        
        Args:
            energies: Energy array
            position: Peak position
            width: Peak width (FWHM)
            amplitude: Peak amplitude
            
        Returns:
            Gaussian profile
        """
        sigma = width / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to sigma
        return amplitude * np.exp(-0.5 * ((energies - position) / sigma) ** 2)
        
    def calculate_spectrum(self, energies: np.ndarray, params: np.ndarray, 
                          theta: float, phi: float) -> np.ndarray:
        """
        Calculate complete spectrum for given measurement geometry.
        
        Args:
            energies: Energy array
            params: Parameter array
            theta: Measurement polar angle in degrees
            phi: Measurement azimuthal angle in degrees
            
        Returns:
            Calculated spectrum
        """
        # Extract global parameters
        i0 = params[0]  # Overall intensity scale
        alpha = np.radians(params[1])  # Sample tilt in radians
        phi_sample = np.radians(params[2])  # Sample azimuth in radians
        
        # Convert measurement angles to radians
        theta_rad = np.radians(theta)
        phi_rad = np.radians(phi)
        
        spectrum = np.zeros_like(energies)
        
        # Loop over transitions
        for i in range(self.n_transitions):
            param_start = self.n_global_params + i * self.n_params_per_transition
            
            # Extract transition parameters
            position = params[param_start]
            width = abs(params[param_start + 1])
            tensor_elements = params[param_start + 2:param_start + 8]
            
            # Create molecular tensor
            molecular_tensor = self.create_molecular_tensor(tensor_elements)
            
            # Transform to film tensor
            film_tensor = self.transform_molecular_tensor(molecular_tensor, alpha)
            
            # Calculate mass absorption
            mass_absorption = self.calculate_mass_absorption(film_tensor, theta_rad, phi_rad)
            
            # Add Gaussian peak with tensor-derived amplitude
            peak_amplitude = i0 * mass_absorption
            spectrum += self.gaussian_peak(energies, position, width, peak_amplitude)
            
        return spectrum
        
    def unpack_parameters(self, params: np.ndarray) -> Dict:
        """
        Unpack parameter array into structured dictionary.
        
        Args:
            params: Flattened parameter array
            
        Returns:
            Dictionary with structured parameters
        """
        result = {
            'i0': params[0],
            'alpha': params[1],
            'phi': params[2],
            'transitions': []
        }
        
        for i in range(self.n_transitions):
            param_start = self.n_global_params + i * self.n_params_per_transition
            
            transition = {
                'position': params[param_start],
                'width': params[param_start + 1],
                'tensor_elements': params[param_start + 2:param_start + 8]
            }
            result['transitions'].append(transition)
            
        return result
        
    def objective_function(self, params: np.ndarray) -> float:
        """
        Objective function for least squares fitting.
        
        Args:
            params: Parameter array
            
        Returns:
            Sum of squared residuals
        """
        total_residual = 0.0
        
        for angle in self.angles:
            # Calculate theoretical spectrum
            calc_spectrum = self.calculate_spectrum(
                self.exp_energies[angle], 
                params, 
                theta=angle,  # Assuming theta = measurement angle
                phi=0  # Assuming phi = 0 for linear polarization
            )
            
            # Calculate residuals
            residuals = self.exp_intensities[angle] - calc_spectrum
            total_residual += np.sum(residuals ** 2)
            
        return total_residual
        
    def setup_bounds(self) -> Bounds:
        """Set up parameter bounds for optimization."""
        lower_bounds = []
        upper_bounds = []
        
        # Global parameter bounds
        lower_bounds.extend([0.001, -90, -180])  # i0, alpha (deg), phi (deg)
        upper_bounds.extend([np.inf, 90, 180])
        
        # Transition parameter bounds
        for i in range(self.n_transitions):
            # Position bounds (±2 eV from initial)
            pos_lower = self.initial_positions[i] - 2.0
            pos_upper = self.initial_positions[i] + 2.0
            
            # Width bounds (0.1 to 5 eV)
            width_lower = 0.1
            width_upper = 5.0
            
            # Tensor element bounds (can be negative)
            tensor_lower = [-10] * 6
            tensor_upper = [10] * 6
            
            lower_bounds.extend([pos_lower, width_lower] + tensor_lower)
            upper_bounds.extend([pos_upper, width_upper] + tensor_upper)
            
        return Bounds(lower_bounds, upper_bounds)
        
    def get_initial_parameters(self) -> np.ndarray:
        """Generate initial parameter array."""
        params = []
        
        # Global parameters: i0, alpha (deg), phi (deg)
        params.extend([1.0, 0.0, 0.0])
        
        # Transition parameters
        for i in range(self.n_transitions):
            # Position and width
            params.extend([self.initial_positions[i], self.initial_widths[i]])
            
            # Tensor elements
            params.extend(self.initial_tensor_elements[i])
            
        return np.array(params)
        
    def fit_spectra(self, method='L-BFGS-B', max_iterations=1000) -> Dict:
        """
        Perform the tensor-based fitting of all spectra simultaneously.
        
        Args:
            method: Optimization method
            max_iterations: Maximum number of iterations
            
        Returns:
            Dictionary with fitting results
        """
        print("Starting tensor-based simultaneous fitting...")
        
        initial_params = self.get_initial_parameters()
        bounds = self.setup_bounds()
        
        print(f"Optimizing {len(initial_params)} parameters...")
        print(f"Initial objective: {self.objective_function(initial_params):.6f}")
        
        # Perform optimization
        result = minimize(
            fun=self.objective_function,
            x0=initial_params,
            method=method,
            bounds=bounds,
            options={'maxiter': max_iterations, 'ftol': 1e-9}
        )
        
        print(f"Final objective: {result.fun:.6f}")
        print(f"Optimization success: {result.success}")
        print(f"Number of iterations: {result.nit}")
        
        # Store results
        self.final_parameters = result.x
        self.fit_result = result
        
        return {
            'result': result,
            'final_params': result.x,
            'objective_value': result.fun,
            'success': result.success,
            'structured_params': self.unpack_parameters(result.x)
        }
        
    def calculate_fit_statistics(self) -> Dict:
        """Calculate detailed fitting statistics."""
        if not hasattr(self, 'final_parameters'):
            raise ValueError("Must run fitting first!")
            
        stats = {}
        total_points = 0
        total_ss_res = 0
        total_ss_tot = 0
        
        for angle in self.angles:
            # Calculate fitted spectrum
            fitted = self.calculate_spectrum(
                self.exp_energies[angle],
                self.final_parameters,
                theta=angle,
                phi=0
            )
            
            # Calculate statistics
            residuals = self.exp_intensities[angle] - fitted
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((self.exp_intensities[angle] - np.mean(self.exp_intensities[angle])) ** 2)
            
            n_points = len(self.exp_intensities[angle])
            rmse = np.sqrt(ss_res / n_points)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            stats[f'angle_{angle}'] = {
                'n_points': n_points,
                'rmse': rmse,
                'r_squared': r_squared,
                'ss_res': ss_res,
                'ss_tot': ss_tot,
                'fitted_spectrum': fitted
            }
            
            total_points += n_points
            total_ss_res += ss_res
            total_ss_tot += ss_tot
            
        # Overall statistics
        stats['overall'] = {
            'total_points': total_points,
            'total_rmse': np.sqrt(total_ss_res / total_points),
            'total_r_squared': 1 - (total_ss_res / total_ss_tot) if total_ss_tot > 0 else 0,
            'n_parameters': len(self.final_parameters)
        }
        
        return stats
        
    def plot_results(self, figsize=(15, 10)):
        """
        Plot experimental vs fitted spectra for all angles.
        
        Args:
            figsize: Figure size tuple
        """
        if not hasattr(self, 'final_parameters'):
            raise ValueError("Must run fitting first!")
            
        n_spectra = len(self.angles)
        n_cols = min(3, n_spectra)
        n_rows = (n_spectra + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n_spectra == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
            
        for i, angle in enumerate(self.angles):
            ax = axes[i] if n_spectra > 1 else axes[0]
            
            # Plot experimental data
            energies = self.exp_energies[angle]
            intensities = self.exp_intensities[angle]
            ax.plot(energies, intensities, 'o-', alpha=0.7, label=f'Exp {angle}°')
            
            # Plot fitted spectrum
            fitted = self.calculate_spectrum(energies, self.final_parameters, 
                                           theta=angle, phi=0)
            ax.plot(energies, fitted, '-', linewidth=2, label=f'Fit {angle}°')
            
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('Intensity')
            ax.set_title(f'Spectrum at θ = {angle}°')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
        # Hide unused subplots
        for j in range(n_spectra, len(axes)):
            axes[j].set_visible(False)
            
        plt.tight_layout()
        plt.show()
        
    def plot_tensor_analysis(self, figsize=(12, 8)):
        """
        Plot analysis of fitted tensor parameters.
        
        Args:
            figsize: Figure size tuple
        """
        if not hasattr(self, 'final_parameters'):
            raise ValueError("Must run fitting first!")
            
        structured_params = self.unpack_parameters(self.final_parameters)
        
        # Extract tensor data
        positions = []
        in_plane_components = []  # Average of Txx and Tyy
        out_plane_components = []  # Tzz
        anisotropies = []
        
        for transition in structured_params['transitions']:
            positions.append(transition['position'])
            tensor = transition['tensor_elements']
            
            in_plane = (tensor[0] + tensor[3]) / 2  # (Txx + Tyy) / 2
            out_plane = tensor[5]  # Tzz
            
            in_plane_components.append(in_plane)
            out_plane_components.append(out_plane)
            
            # Calculate anisotropy ratio
            total = abs(in_plane) + abs(out_plane)
            anisotropy = abs(out_plane) / total if total > 0 else 0
            anisotropies.append(anisotropy)
            
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        
        # Plot tensor components vs energy
        ax1.scatter(positions, in_plane_components, alpha=0.7, label='In-plane (Txx,Tyy)')
        ax1.scatter(positions, out_plane_components, alpha=0.7, label='Out-of-plane (Tzz)')
        ax1.set_xlabel('Energy (eV)')
        ax1.set_ylabel('Tensor Component')
        ax1.set_title('Tensor Components vs Energy')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot anisotropy vs energy
        ax2.scatter(positions, anisotropies, alpha=0.7, color='red')
        ax2.set_xlabel('Energy (eV)')
        ax2.set_ylabel('Out-of-plane Fraction')
        ax2.set_title('Anisotropy vs Energy')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, 1)
        
        # Plot tensor component correlation
        ax3.scatter(in_plane_components, out_plane_components, alpha=0.7)
        ax3.set_xlabel('In-plane Component')
        ax3.set_ylabel('Out-of-plane Component')
        ax3.set_title('Tensor Component Correlation')
        ax3.grid(True, alpha=0.3)
        
        # Plot global parameters
        global_params = ['i0', 'alpha (deg)', 'phi (deg)']
        global_values = [structured_params['i0'], structured_params['alpha'], structured_params['phi']]
        ax4.bar(global_params, global_values)
        ax4.set_title('Global Parameters')
        ax4.set_ylabel('Value')
        
        plt.tight_layout()
        plt.show()
        
    def export_results(self, filename: str):
        """
        Export fitting results to CSV file.
        
        Args:
            filename: Output filename
        """
        if not hasattr(self, 'final_parameters'):
            raise ValueError("Must run fitting first!")
            
        structured_params = self.unpack_parameters(self.final_parameters)
        
        # Create results DataFrame
        results_data = []
        
        for i, transition in enumerate(structured_params['transitions']):
            row = {
                'transition': i,
                'energy': transition['position'],
                'width': transition['width'],
                'Txx': transition['tensor_elements'][0],
                'Txy': transition['tensor_elements'][1],
                'Txz': transition['tensor_elements'][2],
                'Tyy': transition['tensor_elements'][3],
                'Tyz': transition['tensor_elements'][4],
                'Tzz': transition['tensor_elements'][5],
            }
            results_data.append(row)
            
        results_df = pd.DataFrame(results_data)
        
        # Add global parameters as metadata
        metadata_rows = []
        metadata_rows.append(['Global_i0', structured_params['i0']] + [''] * (len(results_df.columns) - 2))
        metadata_rows.append(['Global_alpha_deg', structured_params['alpha']] + [''] * (len(results_df.columns) - 2))
        metadata_rows.append(['Global_phi_deg', structured_params['phi']] + [''] * (len(results_df.columns) - 2))
        metadata_rows.append([''] * len(results_df.columns))  # Empty row
        
        # Create metadata DataFrame
        metadata_df = pd.DataFrame(metadata_rows, columns=results_df.columns)
        
        # Combine and save
        final_df = pd.concat([metadata_df, results_df], ignore_index=True)
        final_df.to_csv(filename, index=False)
        
        print(f"Results exported to {filename}")
        
    def print_fit_summary(self):
        """Print a summary of the fitting results."""
        if not hasattr(self, 'final_parameters'):
            print("No fitting results available. Run fit_spectra() first.")
            return
            
        structured_params = self.unpack_parameters(self.final_parameters)
        stats = self.calculate_fit_statistics()
        
        print("\n" + "="*60)
        print("NEXAFS TENSOR FITTING RESULTS SUMMARY")
        print("="*60)
        
        print(f"\nGlobal Parameters:")
        print(f"  Intensity scale (i0): {structured_params['i0']:.4f}")
        print(f"  Sample tilt (α): {structured_params['alpha']:.2f}°")
        print(f"  Sample azimuth (φ): {structured_params['phi']:.2f}°")
        
        print(f"\nFit Quality:")
        print(f"  Overall R²: {stats['overall']['total_r_squared']:.4f}")
        print(f"  Overall RMSE: {stats['overall']['total_rmse']:.6f}")
        print(f"  Total data points: {stats['overall']['total_points']}")
        
        print(f"\nPer-Spectrum Statistics:")
        for angle in self.angles:
            angle_stats = stats[f'angle_{angle}']
            print(f"  θ = {angle}°: R² = {angle_stats['r_squared']:.4f}, "
                  f"RMSE = {angle_stats['rmse']:.6f}")
        
        print(f"\nTransition Summary:")
        for i, transition in enumerate(structured_params['transitions']):
            tensor = transition['tensor_elements']
            in_plane = (tensor[0] + tensor[3]) / 2
            out_plane = tensor[5]
            total_strength = abs(in_plane) + abs(out_plane)
            anisotropy = abs(out_plane) / total_strength if total_strength > 0 else 0
            
            print(f"  Peak {i+1}: E = {transition['position']:.2f} eV, "
                  f"Width = {transition['width']:.2f} eV, "
                  f"Anisotropy = {anisotropy:.3f}")
        
        print("="*60)


# Example usage and helper functions
def load_and_fit_nexafs_data(param_file: str, exp_file: str, 
                            plot_results: bool = True) -> NEXAFSTensorFitter:
    """
    Convenience function to load data and perform fitting.
    
    Args:
        param_file: Parameter CSV file
        exp_file: Experimental data CSV file
        plot_results: Whether to plot results
        
    Returns:
        Fitted NEXAFSTensorFitter object
    """
    # Initialize fitter
    fitter = NEXAFSTensorFitter(param_file, exp_file)
    
    # Perform fitting
    results = fitter.fit_spectra()
    
    if results['success']:
        print("Fitting completed successfully!")
        fitter.print_fit_summary()
        
        if plot_results:
            fitter.plot_results()
            fitter.plot_tensor_analysis()
    else:
        print("Fitting failed. Check initial parameters and bounds.")
        
    return fitter


# Additional utility functions for data preparation
def create_example_parameter_file(filename: str, n_peaks: int = 5):
    """Create an example parameter file for testing."""
    np.random.seed(42)
    
    data = {
        'E': np.linspace(280, 290, n_peaks),
        'width': np.random.uniform(0.5, 1.5, n_peaks),
        'OS': np.random.uniform(0.1, 1.0, n_peaks),
        'Txx': np.random.uniform(0.3, 0.8, n_peaks),
        'Txy': np.zeros(n_peaks),
        'Txz': np.zeros(n_peaks),
        'Tyy': np.random.uniform(0.3, 0.8, n_peaks),
        'Tyz': np.zeros(n_peaks),
        'Tzz': np.random.uniform(0.1, 0.5, n_peaks)
    }
    
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    print(f"Example parameter file created: {filename}")


def create_synthetic_nexafs_data(param_file: str, output_file: str, 
                                angles: List[float] = None, 
                                energy_range: Tuple[float, float] = (275, 295),
                                n_points: int = 500, noise_level: float = 0.02):
    """
    Create synthetic NEXAFS data for testing the fitting routine.
    
    Args:
        param_file: Parameter file to use for generation
        output_file: Output CSV file for synthetic data
        angles: List of measurement angles
        energy_range: Energy range tuple (min, max)
        n_points: Number of energy points
        noise_level: Relative noise level
    """
    if angles is None:
        angles = [20, 30, 45, 60, 70, 90]
        
    # Load parameters
    params_df = pd.read_csv(param_file)
    
    # Create energy grid
    energies = np.linspace(energy_range[0], energy_range[1], n_points)
    
    # Initialize data storage
    data = {}
    
    # Create a temporary fitter to use its methods
    temp_fitter = NEXAFSTensorFitter.__new__(NEXAFSTensorFitter)
    temp_fitter.params_df = params_df
    temp_fitter.n_transitions = len(params_df)
    temp_fitter.initial_positions = params_df['E'].values
    temp_fitter.initial_widths = params_df['width'].values
    temp_fitter.initial_amplitudes = params_df['OS'].values
    
    # Check if tensor elements exist
    tensor_cols = ['Txx', 'Txy', 'Txz', 'Tyy', 'Tyz', 'Tzz']
    if all(col in params_df.columns for col in tensor_cols):
        temp_fitter.initial_tensor_elements = params_df[tensor_cols].values
    else:
        temp_fitter.initial_tensor_elements = temp_fitter._create_default_tensors()
    
    temp_fitter.setup_tensor_parameters()
    
    # True parameters for generation
    true_params = temp_fitter.get_initial_parameters()
    true_params[0] = 1.0  # i0
    true_params[1] = 15.0  # alpha (degrees)
    true_params[2] = 0.0   # phi (degrees)
    
    # Generate spectra for each angle
    for angle in angles:
        spectrum = temp_fitter.calculate_spectrum(energies, true_params, theta=angle, phi=0)
        
        # Add noise
        noise = np.random.normal(0, noise_level * np.max(spectrum), len(spectrum))
        spectrum_noisy = spectrum + noise
        
        # Store data
        data[f'E_{angle}'] = energies
        data[f'CuPc_CuI_{angle}'] = spectrum_noisy
        
    # Create DataFrame and save
    max_len = max(len(v) for v in data.values())
    for key in data:
        if len(data[key]) < max_len:
            data[key] = np.pad(data[key], (0, max_len - len(data[key])), 'constant', constant_values=np.nan)
            
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    
    print(f"Synthetic NEXAFS data created: {output_file}")
    print(f"True parameters: i0={true_params[0]:.2f}, α={true_params[1]:.1f}°, φ={true_params[2]:.1f}°")


def analyze_tensor_evolution(fitter: NEXAFSTensorFitter, energy_range: Tuple[float, float] = None):
    """
    Analyze how tensor properties evolve with energy.
    
    Args:
        fitter: Fitted NEXAFSTensorFitter object
        energy_range: Energy range for analysis
    """
    if not hasattr(fitter, 'final_parameters'):
        raise ValueError("Must run fitting first!")
        
    structured_params = fitter.unpack_parameters(fitter.final_parameters)
    
    # Extract data
    energies = []
    dichroic_ratios = []
    total_oscillator_strengths = []
    tilt_angles = []  # Molecular tilt relative to surface
    
    for transition in structured_params['transitions']:
        energy = transition['position']
        tensor = transition['tensor_elements']
        
        # Calculate total oscillator strength
        total_os = np.sqrt(np.sum(np.array(tensor)**2))
        
        # Calculate dichroic ratio (out-of-plane / in-plane)
        in_plane = (abs(tensor[0]) + abs(tensor[3])) / 2  # Average |Txx| and |Tyy|
        out_plane = abs(tensor[5])  # |Tzz|
        
        dichroic_ratio = out_plane / in_plane if in_plane > 0 else np.inf
        
        # Estimate molecular tilt angle
        # For a uniaxial system, the ratio Tzz/(Txx+Tyy) relates to molecular orientation
        if (abs(tensor[0]) + abs(tensor[3])) > 0:
            tilt_angle = np.degrees(np.arctan2(out_plane, in_plane))
        else:
            tilt_angle = 90.0
            
        energies.append(energy)
        dichroic_ratios.append(dichroic_ratio)
        total_oscillator_strengths.append(total_os)
        tilt_angles.append(tilt_angle)
        
    # Create comprehensive analysis plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # Dichroic ratio vs energy
    ax1.scatter(energies, dichroic_ratios, s=50, alpha=0.7, c=total_oscillator_strengths, cmap='viridis')
    ax1.set_xlabel('Energy (eV)')
    ax1.set_ylabel('Dichroic Ratio (Out/In-plane)')
    ax1.set_title('Dichroic Ratio vs Energy')
    ax1.grid(True, alpha=0.3)
    cbar1 = plt.colorbar(ax1.collections[0], ax=ax1)
    cbar1.set_label('Total OS')
    
    # Total oscillator strength vs energy
    ax2.scatter(energies, total_oscillator_strengths, s=50, alpha=0.7, c=dichroic_ratios, cmap='plasma')
    ax2.set_xlabel('Energy (eV)')
    ax2.set_ylabel('Total Oscillator Strength')
    ax2.set_title('Oscillator Strength vs Energy')
    ax2.grid(True, alpha=0.3)
    cbar2 = plt.colorbar(ax2.collections[0], ax=ax2)
    cbar2.set_label('Dichroic Ratio')
    
    # Molecular tilt angle vs energy
    ax3.scatter(energies, tilt_angles, s=50, alpha=0.7, c=total_oscillator_strengths, cmap='cool')
    ax3.set_xlabel('Energy (eV)')
    ax3.set_ylabel('Effective Tilt Angle (degrees)')
    ax3.set_title('Molecular Orientation vs Energy')
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 90)
    cbar3 = plt.colorbar(ax3.collections[0], ax=ax3)
    cbar3.set_label('Total OS')
    
    # Summary statistics
    ax4.hist(dichroic_ratios, bins=10, alpha=0.7, label='Dichroic Ratios')
    ax4.axvline(np.mean(dichroic_ratios), color='red', linestyle='--', label=f'Mean: {np.mean(dichroic_ratios):.2f}')
    ax4.set_xlabel('Dichroic Ratio')
    ax4.set_ylabel('Count')
    ax4.set_title('Distribution of Dichroic Ratios')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Print summary statistics
    print("\n" + "="*50)
    print("TENSOR EVOLUTION ANALYSIS")
    print("="*50)
    print(f"Energy range: {min(energies):.2f} - {max(energies):.2f} eV")
    print(f"Mean dichroic ratio: {np.mean(dichroic_ratios):.3f} ± {np.std(dichroic_ratios):.3f}")
    print(f"Mean total OS: {np.mean(total_oscillator_strengths):.3f} ± {np.std(total_oscillator_strengths):.3f}")
    print(f"Mean tilt angle: {np.mean(tilt_angles):.1f}° ± {np.std(tilt_angles):.1f}°")
    print(f"Sample tilt (α): {structured_params['alpha']:.1f}°")
    print("="*50)


def compare_fitting_methods(param_file: str, exp_file: str) -> Dict:
    """
    Compare tensor-based fitting with simple Gaussian fitting.
    
    Args:
        param_file: Parameter file
        exp_file: Experimental data file
        
    Returns:
        Dictionary with comparison results
    """
    print("Comparing tensor-based vs simple Gaussian fitting...")
    
    # Tensor-based fitting
    print("\n--- Tensor-based fitting ---")
    tensor_fitter = NEXAFSTensorFitter(param_file, exp_file)
    tensor_results = tensor_fitter.fit_spectra()
    tensor_stats = tensor_fitter.calculate_fit_statistics()
    
    # Simple Gaussian fitting (using original class)
    print("\n--- Simple Gaussian fitting ---")
    try:
        simple_fitter = MultiSpectrumGaussianFitter(param_file, exp_file)
        simple_results = simple_fitter.fit_all_stages()
        simple_stats = simple_fitter.calculate_fit_statistics()
    except:
        print("Simple Gaussian fitter not available for comparison")
        simple_results = None
        simple_stats = None
    
    # Compare results
    comparison = {
        'tensor_based': {
            'r_squared': tensor_stats['overall']['total_r_squared'],
            'rmse': tensor_stats['overall']['total_rmse'],
            'n_parameters': len(tensor_fitter.final_parameters),
            'success': tensor_results['success']
        }
    }
    
    if simple_results and simple_stats:
        comparison['simple_gaussian'] = {
            'r_squared': simple_stats['overall']['total_r_squared'],
            'rmse': simple_stats['overall']['total_rmse'],
            'n_parameters': len(simple_fitter.final_parameters),
            'success': simple_results['stage3']['success']
        }
        
        print(f"\n{'Method':<20} {'R²':<10} {'RMSE':<12} {'N_params':<12}")
        print("-" * 54)
        print(f"{'Tensor-based':<20} {comparison['tensor_based']['r_squared']:<10.4f} "
              f"{comparison['tensor_based']['rmse']:<12.6f} {comparison['tensor_based']['n_parameters']:<12}")
        if 'simple_gaussian' in comparison:
            print(f"{'Simple Gaussian':<20} {comparison['simple_gaussian']['r_squared']:<10.4f} "
                  f"{comparison['simple_gaussian']['rmse']:<12.6f} {comparison['simple_gaussian']['n_parameters']:<12}")
    
    return comparison
