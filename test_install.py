#!/usr/bin/env python3
"""
Test that the package can be installed and imported correctly
"""

def test_imports():
    """Test basic package imports"""
    try:
        import radiosonde_teamxuk
        print(f"✓ Package imported successfully")
        print(f"  Version: {radiosonde_teamxuk.__version__}")
        
        from radiosonde_teamxuk import read_edt_file, save_to_ncas_netcdf
        print(f"✓ process_radiosondes functions imported")
        
        from radiosonde_teamxuk import create_quicklook
        print(f"✓ generate_quicklooks functions imported")
        
        # Check if metadata files are accessible
        try:
            from importlib.resources import files
            pkg_files = files('radiosonde_teamxuk')
            
            metadata_files = [
                'metadata_ncas-radiosonde-1.json',
                'metadata_ncas-radiosonde-2.json',
                'AMF_product_sonde_variable.json'
            ]
            
            for fname in metadata_files:
                if (pkg_files / fname).is_file():
                    print(f"✓ Found metadata file: {fname}")
                else:
                    print(f"⚠ Missing metadata file: {fname}")
        except Exception as e:
            print(f"⚠ Could not check metadata files: {e}")
        
        print("\n✓ All imports successful!")
        return True
        
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False

if __name__ == "__main__":
    success = test_imports()
    exit(0 if success else 1)
