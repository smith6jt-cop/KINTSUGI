#!/usr/bin/env python3
"""
Simple Environment Checker for KINTSUGI Registration Workflow
============================================================

This script checks if your environment is ready for the registration workflow
and provides setup instructions.

Author: KINTSUGI Image Processing Pipeline
"""

import os
import sys
from pathlib import Path

def check_basic_environment():
    """Check basic Python environment."""
    print("=== Basic Environment Check ===")
    print(f"Python version: {sys.version}")
    print(f"Python executable: {sys.executable}")
    
    # Check if we're in a conda environment
    if 'CONDA_DEFAULT_ENV' in os.environ:
        print(f"Conda environment: {os.environ['CONDA_DEFAULT_ENV']}")
    else:
        print("Not in a conda environment")
    
    # Check current directory
    current_dir = os.getcwd()
    print(f"Current directory: {current_dir}")
    
    # Check if we're in the KINTSUGI directory
    if 'KINTSUGI' in current_dir:
        print("✓ Running from KINTSUGI directory")
    else:
        print("⚠️  Not in KINTSUGI directory")
    
    return True

def check_kintsugi_structure():
    """Check KINTSUGI directory structure."""
    print("\n=== KINTSUGI Structure Check ===")
    
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    required_files = [
        'env.yml',
        'README.md', 
        'notebooks/Kreg/__init__.py',
        'notebooks/Kreg/registration.py',
        'notebooks/Kreg/slide_io.py'
    ]
    
    missing_files = []
    for file_path in required_files:
        full_path = os.path.join(base_dir, file_path)
        if os.path.exists(full_path):
            print(f"✓ {file_path}")
        else:
            print(f"✗ {file_path}")
            missing_files.append(file_path)
    
    if missing_files:
        print(f"\n⚠️  Missing {len(missing_files)} required files")
        return False
    else:
        print("✓ All required KINTSUGI files found")
        return True

def check_workflow_files():
    """Check if registration workflow files exist."""
    print("\n=== Registration Workflow Files ===")
    
    workflow_files = [
        'Image_Registration_Workflow.ipynb',
        'registration_workflow.py',
        'config_example.json',
        'sample_channel_names.txt',
        'README_Registration.md'
    ]
    
    notebooks_dir = os.path.dirname(os.path.abspath(__file__))
    
    for file_name in workflow_files:
        file_path = os.path.join(notebooks_dir, file_name)
        if os.path.exists(file_path):
            size = os.path.getsize(file_path)
            print(f"✓ {file_name} ({size} bytes)")
        else:
            print(f"✗ {file_name}")

def print_setup_instructions():
    """Print setup instructions."""
    print("\n" + "="*60)
    print("SETUP INSTRUCTIONS")
    print("="*60)
    
    print("\n1. ENVIRONMENT SETUP:")
    print("   The KINTSUGI registration workflow requires the full KINTSUGI")
    print("   conda environment with all scientific Python packages.")
    print()
    print("   To set up the environment:")
    print("   a) Open Anaconda Prompt (Windows) or Terminal (Mac/Linux)")
    print("   b) Navigate to the KINTSUGI directory")
    print("   c) Run: conda env create -f env.yml")
    print("   d) Run: conda activate KINTSUGI")
    print("   e) Run: code . (to open VS Code from the activated environment)")
    
    print("\n2. REQUIRED FILES:")
    print("   Place your image files in the 'data/' directory:")
    print("   • For PANCREAS: any TIFF files with 'pancreas' in the name")
    print("   • For THYMUS: ThymusDAPI_Xenium_downsampled2x_EDF.ome.tiff")
    print("                 and THYMUS.ome.tiff")
    print("   • Optional: channel names text files")
    
    print("\n3. RUNNING THE WORKFLOW:")
    print("   Once the environment is set up, you can run:")
    print("   • python registration_workflow.py --help (for command line)")
    print("   • Open Image_Registration_Workflow.ipynb in Jupyter/VS Code")
    
    print("\n4. TROUBLESHOOTING:")
    print("   If you encounter issues:")
    print("   • Ensure you're in the KINTSUGI conda environment")
    print("   • Check that env.yml was used to create the environment")
    print("   • Verify Java and Maven are properly installed")
    print("   • Restart your terminal/VS Code if needed")

def create_directories():
    """Create necessary directories."""
    print("\n=== Creating Directories ===")
    
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    directories = [
        os.path.join(base_dir, 'data'),
        os.path.join(base_dir, 'outputs')
    ]
    
    for directory in directories:
        if not os.path.exists(directory):
            try:
                os.makedirs(directory, exist_ok=True)
                print(f"✓ Created: {directory}")
            except Exception as e:
                print(f"✗ Failed to create {directory}: {e}")
        else:
            print(f"✓ Exists: {directory}")

def main():
    """Main function."""
    print("KINTSUGI Registration Workflow - Environment Checker")
    print("=" * 60)
    
    # Check basic environment
    check_basic_environment()
    
    # Check KINTSUGI structure
    structure_ok = check_kintsugi_structure()
    
    # Check workflow files
    check_workflow_files()
    
    # Create necessary directories
    create_directories()
    
    # Print instructions
    print_setup_instructions()
    
    if structure_ok:
        print("\n✅ KINTSUGI structure looks good!")
        print("   Follow the setup instructions above to prepare your environment.")
    else:
        print("\n❌ KINTSUGI structure issues detected.")
        print("   Please ensure you have the complete KINTSUGI repository.")
    
    print(f"\n📁 Current working directory: {os.getcwd()}")
    print("   Make sure to run the workflow from the notebooks/ directory")
    print("   in the activated KINTSUGI conda environment.")

if __name__ == "__main__":
    main()