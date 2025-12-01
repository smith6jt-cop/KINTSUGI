#!/usr/bin/env python3
"""
Enhanced Registration Workflow with Better Error Handling
=========================================================

This script provides an enhanced version of the image registration workflow
with improved error handling, progress reporting, and troubleshooting capabilities.

Author: KINTSUGI Image Processing Pipeline
"""

import os
import sys
import traceback
import time
from pathlib import Path

def check_environment():
    """Check if the environment is properly set up for registration."""
    print("=== Environment Check ===")
    
    # Check Python version
    print(f"Python version: {sys.version}")
    
    # Check required modules
    required_modules = [
        'numpy', 'tifffile', 'PIL', 'pathlib', 'itertools'
    ]
    
    missing_modules = []
    for module in required_modules:
        try:
            __import__(module)
            print(f"✓ {module} available")
        except ImportError:
            print(f"✗ {module} missing")
            missing_modules.append(module)
    
    # Check KINTSUGI specific modules
    kreg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Kreg')
    if os.path.exists(kreg_path):
        print(f"✓ Kreg module found at: {kreg_path}")
        
        # Try to import Kreg modules
        try:
            sys.path.append(kreg_path)
            from Kreg import registration
            print("✓ Kreg.registration imported successfully")
        except Exception as e:
            print(f"✗ Failed to import Kreg.registration: {e}")
            return False
    else:
        print(f"✗ Kreg module not found at: {kreg_path}")
        return False
    
    if missing_modules:
        print(f"\n✗ Missing modules: {', '.join(missing_modules)}")
        print("Please ensure you're running in the KINTSUGI conda environment:")
        print("  conda activate KINTSUGI")
        return False
    
    print("✓ Environment check passed")
    return True


def create_sample_data():
    """Create sample data files for testing."""
    print("\n=== Creating Sample Data ===")
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    
    # Create sample channel names file
    sample_channels = [
        "DAPI", "CD3", "CD4", "CD8", "CD20", "CD31", "CD68", "Ki67",
        "PanCK", "CD45", "Vimentin", "CD21", "CD34", "CD44", "FoxP3",
        "CD163", "CD107a", "HLADR", "CollagenIV", "ECAD", "CD45RO",
        "CD11c", "CD35", "Lyve1", "SMActin", "Podoplanin", "CD15",
        "CD5", "CD1c", "Cytokeratin"
    ]
    
    # Extend to 62 channels
    while len(sample_channels) < 62:
        sample_channels.append(f"Channel{len(sample_channels)+1}")
    
    # Write sample channel names
    channel_file = os.path.join(data_dir, "sample_channel_names.txt")
    with open(channel_file, 'w') as f:
        for channel in sample_channels:
            f.write(f"{channel}\n")
    
    print(f"✓ Created sample channel names file: {channel_file}")
    
    # Create simple test configuration
    config = {
        "data_dir": data_dir,
        "output_dir": os.path.join(base_dir, "outputs"),
        "channel_names_file": channel_file,
        "max_processed_image_dim_px": 1000,
        "max_non_rigid_registration_dim_px": 2000
    }
    
    return config


def safe_import_dependencies():
    """Safely import all required dependencies."""
    try:
        import numpy as np
        import tifffile
        from PIL import Image
        Image.MAX_IMAGE_PIXELS = 2000000000
        
        # Add Kreg to path and import
        kreg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Kreg')
        if kreg_path not in sys.path:
            sys.path.append(kreg_path)
        
        from Kreg import registration
        from Kreg import slide_io
        
        return {
            'np': np,
            'tifffile': tifffile,
            'Image': Image,
            'registration': registration,
            'slide_io': slide_io
        }
    except Exception as e:
        print(f"Error importing dependencies: {e}")
        return None


def diagnose_issues(config):
    """Diagnose common issues and provide solutions."""
    print("\n=== Diagnostic Report ===")
    
    issues_found = []
    solutions = []
    
    # Check data directory
    if not os.path.exists(config['data_dir']):
        issues_found.append("Data directory does not exist")
        solutions.append(f"Create data directory: mkdir -p {config['data_dir']}")
    
    # Check for image files
    image_extensions = ['.tiff', '.tif', '.ome.tiff', '.ome.tif']
    image_files = []
    
    if os.path.exists(config['data_dir']):
        for file in os.listdir(config['data_dir']):
            if any(file.lower().endswith(ext) for ext in image_extensions):
                image_files.append(file)
    
    if len(image_files) < 2:
        issues_found.append(f"Not enough image files found (need at least 2, found {len(image_files)})")
        solutions.append("Add image files to the data directory")
    
    # Check disk space
    try:
        import shutil
        total, used, free = shutil.disk_usage(config.get('output_dir', '.'))
        free_gb = free / (1024**3)
        if free_gb < 5:
            issues_found.append(f"Low disk space ({free_gb:.1f} GB free)")
            solutions.append("Free up disk space (recommend at least 5 GB)")
    except:
        pass
    
    # Print results
    if issues_found:
        print("Issues found:")
        for i, issue in enumerate(issues_found, 1):
            print(f"  {i}. {issue}")
        
        print("\nSuggested solutions:")
        for i, solution in enumerate(solutions, 1):
            print(f"  {i}. {solution}")
    else:
        print("✓ No obvious issues detected")
    
    return len(issues_found) == 0


def test_registration_pipeline():
    """Test the registration pipeline with minimal data."""
    print("\n=== Testing Registration Pipeline ===")
    
    try:
        # Import dependencies
        deps = safe_import_dependencies()
        if deps is None:
            return False
        
        np = deps['np']
        registration = deps['registration']
        
        print("✓ Dependencies imported successfully")
        
        # Create minimal test images
        test_dir = "test_registration"
        os.makedirs(test_dir, exist_ok=True)
        
        # Create simple test images (100x100 pixels)
        test_img1 = np.random.randint(0, 65535, (100, 100), dtype=np.uint16)
        test_img2 = np.roll(test_img1, 5, axis=0)  # Shift image slightly
        
        deps['tifffile'].imwrite(os.path.join(test_dir, "test1.tiff"), test_img1)
        deps['tifffile'].imwrite(os.path.join(test_dir, "test2.tiff"), test_img2)
        
        print("✓ Test images created")
        
        # Test VALIS initialization
        results_dir = os.path.join(test_dir, "results")
        os.makedirs(results_dir, exist_ok=True)
        
        # This is a minimal test - just check if Valis can be initialized
        try:
            registrar = registration.Valis(
                src_dir=test_dir,
                dst_dir=results_dir,
                max_processed_image_dim_px=50,  # Very small for testing
                max_non_rigid_registration_dim_px=50
            )
            print("✓ VALIS registrar initialized successfully")
            
            # Clean up
            import shutil
            shutil.rmtree(test_dir)
            
            return True
            
        except Exception as e:
            print(f"✗ VALIS initialization failed: {e}")
            
            # Clean up
            import shutil
            shutil.rmtree(test_dir, ignore_errors=True)
            
            return False
            
    except Exception as e:
        print(f"✗ Pipeline test failed: {e}")
        traceback.print_exc()
        return False


def main():
    """Main function for enhanced registration workflow."""
    print("KINTSUGI Enhanced Registration Workflow")
    print("=" * 50)
    
    # Step 1: Check environment
    if not check_environment():
        print("\n❌ Environment check failed. Please fix the issues above.")
        return False
    
    # Step 2: Create sample data
    config = create_sample_data()
    
    # Step 3: Diagnose issues
    if not diagnose_issues(config):
        print("\n⚠️  Issues detected. Please address them before running the full workflow.")
    
    # Step 4: Test pipeline
    if test_registration_pipeline():
        print("\n✅ Registration pipeline test passed!")
        print("\nYou can now run the full registration workflow:")
        print("  1. Place your image files in the data/ directory")
        print("  2. Run: python registration_workflow.py")
        print("  3. Or use the Jupyter notebook: Image_Registration_Workflow.ipynb")
        return True
    else:
        print("\n❌ Registration pipeline test failed.")
        print("\nTroubleshooting steps:")
        print("  1. Ensure you're in the KINTSUGI conda environment")
        print("  2. Check that all required packages are installed")
        print("  3. Verify that Java is properly configured")
        print("  4. Try restarting your Python session")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)