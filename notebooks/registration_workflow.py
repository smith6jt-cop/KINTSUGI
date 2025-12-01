#!/usr/bin/env python3
"""
KINTSUGI Image Registration and Merging Workflow
=================================================

This script provides a workflow for registering and merging PANCREAS and THYMUS images,
ensuring proper cropping, channel naming, and output as pyramidal TIFF files compatible with QuPath.

Usage:
    python registration_workflow.py [--dataset DATASET] [--config CONFIG_FILE]
    
Arguments:
    --dataset: Specify 'pancreas', 'thymus', or 'both' (default: both)
    --config: Path to configuration file (optional)

Author: KINTSUGI Image Processing Pipeline
"""

import os
import sys
import argparse
import json
import numpy as np
import tifffile
from pathlib import Path
import traceback
import shutil
from itertools import chain

# Add the Kreg module to path
kreg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Kreg')
if kreg_path not in sys.path:
    sys.path.append(kreg_path)

try:
    from Kreg import registration
    from Kreg import slide_io
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = 2000000000
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Please ensure the Kreg module is available and all dependencies are installed.")
    sys.exit(1)


def safe_channel_names(path):
    """Safely load channel names from a file.
    
    Parameters
    ----------
    path : str or Path
        Path to the channel names file
        
    Returns
    -------
    list
        List of channel names, or default names if file not found
    """
    try:
        if path and os.path.exists(path):
            with open(path, 'r') as f:
                channel_names = [line.strip() for line in f.readlines() if line.strip()]
            print(f"Loaded {len(channel_names)} channel names from {path}")
            return channel_names
        else:
            print(f"Channel names file not found at {path}, using default names")
            # Generate default 62-channel names
            return [f"C{i+1}" for i in range(62)]
    except Exception as e:
        print(f"Error loading channel names: {e}")
        return [f"C{i+1}" for i in range(62)]


def safe_array(arr):
    """Ensure the array is in the correct shape and format.
    
    Parameters
    ----------
    arr : numpy.ndarray
        Input array
        
    Returns
    -------
    numpy.ndarray
        Array in correct format
    """
    try:
        if arr is None:
            return None
        
        # Ensure array is numpy array
        if not isinstance(arr, np.ndarray):
            arr = np.asarray(arr)
        
        # Ensure proper dtype for image data
        if arr.dtype == np.bool_:
            arr = arr.astype(np.uint8) * 255
        elif np.issubdtype(arr.dtype, np.integer) and arr.dtype != np.uint16:
            # Convert to uint16 for better compatibility
            if arr.max() <= 255:
                arr = arr.astype(np.uint16) * 257  # Scale 8-bit to 16-bit
            else:
                arr = arr.astype(np.uint16)
        
        return arr
    except Exception as e:
        print(f"Error processing array: {e}")
        return arr


def check_file_exists(filepath, description="File"):
    """Check if file exists and provide helpful error message.
    
    Parameters
    ----------
    filepath : str or Path
        Path to check
    description : str
        Description of the file for error messages
        
    Returns
    -------
    bool
        True if file exists, False otherwise
    """
    if filepath and os.path.exists(filepath):
        print(f"✓ {description} found: {filepath}")
        return True
    else:
        print(f"✗ {description} not found: {filepath}")
        return False


def register_and_merge_images(config):
    """Register and merge images according to the KINTSUGI workflow.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary with keys:
        - name: str, name of the dataset
        - fixed_image: str, path to fixed image
        - moving_image: str, path to moving image  
        - channel_names_file: str, path to channel names file
        - output_file: str, path for output file
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        print(f"\n=== Starting {config['name']} Registration and Merging ===")
        
        # Check if files exist
        if not check_file_exists(config['fixed_image'], "Fixed image"):
            return False
        if not check_file_exists(config['moving_image'], "Moving image"):
            return False
            
        # Step 1: Load channel names
        print("\n1. Loading channel names...")
        channel_names = safe_channel_names(config.get('channel_names_file'))
        
        # Step 2: Create temporary directory for registration
        temp_dir = os.path.join(os.path.dirname(config['output_file']), f"temp_{config['name']}_registration")
        os.makedirs(temp_dir, exist_ok=True)
        
        # Copy images to temp directory for registration
        temp_fixed = os.path.join(temp_dir, "fixed.ome.tiff")
        temp_moving = os.path.join(temp_dir, "moving.ome.tiff")
        
        if not os.path.exists(temp_fixed):
            shutil.copy2(config['fixed_image'], temp_fixed)
        if not os.path.exists(temp_moving):
            shutil.copy2(config['moving_image'], temp_moving)
        
        # Step 3: Initialize registration
        print("\n2. Initializing VALIS registration...")
        results_dir = os.path.join(temp_dir, "results")
        os.makedirs(results_dir, exist_ok=True)
        
        registrar = registration.Valis(
            src_dir=temp_dir,
            dst_dir=results_dir,
            max_processed_image_dim_px=1000,  # Adjust based on image size
            max_non_rigid_registration_dim_px=2000,
            crop="overlap",  # Crop to overlapping regions
            compose_non_rigid=True
        )
        
        # Step 4: Perform registration
        print("\n3. Performing registration...")
        rigid_registrar, non_rigid_registrar, error_df = registrar.register()
        
        print(f"Registration completed. Error statistics:")
        if error_df is not None and not error_df.empty:
            print(f"  Mean error: {error_df['error'].mean():.2f}")
            print(f"  Max error: {error_df['error'].max():.2f}")
        
        # Step 5: Warp and merge slides
        print("\n4. Warping and merging slides...")
        
        # Create channel name dictionary
        slide_names = list(registrar.slide_dict.keys())
        channels_per_slide = len(channel_names) // len(slide_names) if len(slide_names) > 0 else len(channel_names)
        
        channel_name_dict = {}
        for i, slide_name in enumerate(slide_names):
            start_idx = i * channels_per_slide
            end_idx = min(start_idx + channels_per_slide, len(channel_names))
            channel_name_dict[slide_name] = channel_names[start_idx:end_idx]
        
        # Perform the merge
        warped_slide = registrar.warp_and_merge_slides(
            dst_f=None,  # Don't save yet, we'll process further
            crop="overlap",
            channel_name_dict=channel_name_dict,
            pyramid=True,
            compression="lzw"
        )
        
        merged_array, all_channel_names, ome_xml = warped_slide
        
        # Step 6: Process merged array to ensure 62 channels
        print("\n5. Processing merged array...")
        
        # Convert pyvips image to numpy array if needed
        if hasattr(merged_array, 'numpy'):
            # For pyvips images
            merged_np = merged_array.numpy()
        else:
            merged_np = np.asarray(merged_array)
        
        merged_np = safe_array(merged_np)
        
        # Ensure we have the right dimensions
        if merged_np.ndim == 3:
            h, w, c = merged_np.shape
            # Transpose to (channels, height, width) for tifffile
            merged_np = np.transpose(merged_np, (2, 0, 1))
        elif merged_np.ndim == 2:
            h, w = merged_np.shape
            c = 1
            merged_np = merged_np[np.newaxis, :, :]
        else:
            c, h, w = merged_np.shape
        
        # Pad or crop to exactly 62 channels
        if c < 62:
            # Pad with zeros
            padding = np.zeros((62 - c, h, w), dtype=merged_np.dtype)
            merged_np = np.concatenate([merged_np, padding], axis=0)
            # Extend channel names
            while len(all_channel_names) < 62:
                all_channel_names.append(f"Empty{len(all_channel_names)+1}")
        elif c > 62:
            # Crop to first 62 channels
            merged_np = merged_np[:62, :, :]
            all_channel_names = all_channel_names[:62]
        
        print(f"Final array shape: {merged_np.shape}")
        print(f"Number of channels: {len(all_channel_names)}")
        
        # Step 7: Create OME-XML metadata
        print("\n6. Creating OME-XML metadata...")
        
        ome_xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"
     xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
     xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">
    <Image ID="Image:0" Name="{config['name']}_registered_62ch">
        <Pixels BigEndian="false"
                DimensionOrder="XYCZT"
                ID="Pixels:0"
                Interleaved="false"
                SignificantBits="16"
                SizeC="62"
                SizeT="1"
                SizeX="{w}"
                SizeY="{h}"
                SizeZ="1"
                Type="uint16">
'''
        
        # Add channel information
        for i, channel_name in enumerate(all_channel_names):
            ome_xml += f'            <Channel ID="Channel:0:{i}" Name="{channel_name}" SamplesPerPixel="1"/>\n'
        
        ome_xml += '''        </Pixels>
    </Image>
</OME>'''
        
        # Step 8: Save pyramidal TIFF
        print("\n7. Saving pyramidal TIFF...")
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(config['output_file']), exist_ok=True)
        
        with tifffile.TiffWriter(config['output_file'], bigtiff=True) as tif:
            tif.write(
                merged_np,
                compression='lzw',
                description=ome_xml,
                metadata={'axes': 'CYX'},
                tile=(512, 512),
                pyramid=True
            )
        
        print(f"\n✓ Successfully saved {config['name']} registered image to:")
        print(f"  {config['output_file']}")
        print(f"  Size: {os.path.getsize(config['output_file']) / (1024**3):.2f} GB")
        print(f"  Dimensions: {merged_np.shape}")
        print(f"  Channels: {len(all_channel_names)}")
        
        # Clean up temporary directory
        try:
            shutil.rmtree(temp_dir)
            print(f"  Cleaned up temporary directory")
        except:
            print(f"  Warning: Could not clean up temporary directory: {temp_dir}")
        
        return True
        
    except Exception as e:
        print(f"\n✗ Error during {config['name']} registration: {e}")
        print("Full traceback:")
        traceback.print_exc()
        return False
    
    finally:
        # Always try to kill JVM if it was started
        try:
            registration.kill_jvm()
        except:
            pass


def load_config(config_file):
    """Load configuration from JSON file."""
    try:
        with open(config_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading config file {config_file}: {e}")
        return None


def get_default_config():
    """Get default configuration."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")
    output_dir = os.path.join(base_dir, "outputs")
    
    return {
        "pancreas": {
            "name": "PANCREAS",
            "fixed_image": None,  # To be auto-detected
            "moving_image": None,  # To be auto-detected
            "channel_names_file": None,  # To be auto-detected
            "output_file": os.path.join(output_dir, "PANCREAS_registered_62ch_ImageJ.tiff")
        },
        "thymus": {
            "name": "THYMUS",
            "fixed_image": os.path.join(data_dir, "ThymusDAPI_Xenium_downsampled2x_EDF.ome.tiff"),
            "moving_image": os.path.join(data_dir, "THYMUS.ome.tiff"),
            "channel_names_file": os.path.join(data_dir, "thymus_channel_names.txt"),
            "output_file": os.path.join(output_dir, "THYMUS_registered_62ch_ImageJ.tiff")
        }
    }


def auto_detect_pancreas_files(data_dir):
    """Auto-detect PANCREAS files in data directory."""
    pancreas_files = []
    if os.path.exists(data_dir):
        for file in os.listdir(data_dir):
            if 'pancreas' in file.lower() and file.endswith(('.tiff', '.tif', '.ome.tiff')):
                pancreas_files.append(os.path.join(data_dir, file))
    
    return pancreas_files


def main():
    parser = argparse.ArgumentParser(description='KINTSUGI Image Registration and Merging Workflow')
    parser.add_argument('--dataset', choices=['pancreas', 'thymus', 'both'], default='both',
                        help='Which dataset to process (default: both)')
    parser.add_argument('--config', type=str, help='Path to configuration JSON file')
    parser.add_argument('--data-dir', type=str, help='Override data directory path')
    parser.add_argument('--output-dir', type=str, help='Override output directory path')
    
    args = parser.parse_args()
    
    # Load configuration
    if args.config:
        config = load_config(args.config)
        if config is None:
            sys.exit(1)
    else:
        config = get_default_config()
    
    # Override directories if specified
    if args.data_dir:
        data_dir = args.data_dir
    else:
        data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    
    if args.output_dir:
        output_dir = args.output_dir
        # Update output paths
        config['pancreas']['output_file'] = os.path.join(output_dir, "PANCREAS_registered_62ch_ImageJ.tiff")
        config['thymus']['output_file'] = os.path.join(output_dir, "THYMUS_registered_62ch_ImageJ.tiff")
    else:
        output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    print("=== KINTSUGI Image Registration and Merging Workflow ===")
    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Processing: {args.dataset}")
    
    success_count = 0
    total_count = 0
    
    # Process PANCREAS
    if args.dataset in ['pancreas', 'both']:
        total_count += 1
        print("\n" + "="*60)
        print("PANCREAS Registration Workflow")
        print("="*60)
        
        # Auto-detect PANCREAS files if not specified
        if not config['pancreas']['fixed_image']:
            pancreas_files = auto_detect_pancreas_files(data_dir)
            print(f"Found {len(pancreas_files)} PANCREAS files in data directory:")
            for f in pancreas_files:
                print(f"  - {f}")
            
            if len(pancreas_files) >= 2:
                config['pancreas']['fixed_image'] = pancreas_files[0]
                config['pancreas']['moving_image'] = pancreas_files[1]
                
                # Look for channel names file
                channel_files = [f for f in os.listdir(data_dir) if 'pancreas' in f.lower() and f.endswith('.txt')]
                if channel_files:
                    config['pancreas']['channel_names_file'] = os.path.join(data_dir, channel_files[0])
        
        if config['pancreas']['fixed_image'] and config['pancreas']['moving_image']:
            if register_and_merge_images(config['pancreas']):
                success_count += 1
                print("✓ PANCREAS registration completed successfully!")
            else:
                print("✗ PANCREAS registration failed.")
        else:
            print("✗ Not enough PANCREAS files found. Need at least 2 image files.")
            print("  Please place PANCREAS image files in the data directory.")
    
    # Process THYMUS
    if args.dataset in ['thymus', 'both']:
        total_count += 1
        print("\n" + "="*60)
        print("THYMUS Registration Workflow")
        print("="*60)
        
        if register_and_merge_images(config['thymus']):
            success_count += 1
            print("✓ THYMUS registration completed successfully!")
        else:
            print("✗ THYMUS registration failed.")
            print("  Please ensure the following files exist in the data directory:")
            print(f"  - {config['thymus']['fixed_image']}")
            print(f"  - {config['thymus']['moving_image']}")
            print(f"  - {config['thymus']['channel_names_file']} (optional)")
    
    # Summary
    print("\n" + "="*60)
    print("WORKFLOW SUMMARY")
    print("="*60)
    print(f"Successfully processed: {success_count}/{total_count} datasets")
    
    if success_count > 0:
        print("\nOutput files:")
        for dataset_name in ['pancreas', 'thymus']:
            if args.dataset in [dataset_name, 'both']:
                output_file = config[dataset_name]['output_file']
                if os.path.exists(output_file):
                    file_size = os.path.getsize(output_file) / (1024**3)
                    print(f"  ✓ {os.path.basename(output_file)} ({file_size:.2f} GB)")
                    print(f"    {output_file}")
        
        print("\nThe registered and merged images are now compatible with QuPath.")
        print("Each file contains 62 channels in pyramidal TIFF format with OME metadata.")
    
    return success_count == total_count


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)