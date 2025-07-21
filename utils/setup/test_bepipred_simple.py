#!/usr/bin/env python3
"""
Simple test for BepiPred 3.0 - direct execution without environment activation
"""

import subprocess
import tempfile
import os
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def test_bepipred_direct():
    """Test BepiPred by directly calling the Python script"""
    
    # Test sequence (dnaK from your data)
    test_sequence = "MSKGPAVGIDLGTTYSCVGVFQHGKVEIIANDQGNRTTPSYVAFTDTERLIGDAAKNQVAMNPTNTVFDAKRLIGRRFDDAVVQSDMKHWPFMVVNDAGRPKVQVEYKGETKSFYPEEVSSMVLTKMKEIAEAYLGKTVTNAVVTVPAYFNDSQRQATKDAGTIAGLNVLRIINEPTAAAIAYGLDKKVGAERNVLIFDLGGGTFDVSILTIEDGIFEVKSTAGDTHLGGEDFDNRMVNHFIAEFKRKHKKDISENKRAVRRLRTACERAKRTLSSSTQASIEIDSLYEGIDFYTSITRARFEELNADLFRGTLDPVEKALRDAKLDKSQIHDIVLVGGSTRIPKIQKLLQDFFNGKELNKSINPDEAVAYGAAVQAAILSGDKSENVQDLLLLDVAPLSLGLETAGGVMTVLIKRNTTIPTKQTQIFTTYSDNQPGVLIQVYEGERAMTKDNNLLGRFELSGIPPAPRGVPQIEVTFDIDANGILNVTATDKSTGKANKITITNDKGRLSKEEIERMVQEAEKYKAEDEVQRERVSAKNALESYAFNMKSAVEDEGLKGKISEADKKKVLDKCQEVISWLDANTLAEKDEFEHKRKELEQVCNPIISGLYQGAGGPGPGGFGAQGPKGGSGSGPTIEEVD"
    
    bepipred_path = Path("tools/BepiPred3.0")
    bepipred_script = bepipred_path / "bepipred3_CLI.py"
    
    print(f"Testing BepiPred 3.0...")
    print(f"BepiPred script: {bepipred_script}")
    print(f"Sequence length: {len(test_sequence)}")
    print()
    
    if not bepipred_script.exists():
        print(f"❌ BepiPred script not found: {bepipred_script}")
        return False
    
    try:
        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">test_dnaK\n{test_sequence}\n")
            fasta_file = f.name
        
        output_dir = "test_bepipred_output"
        os.makedirs(output_dir, exist_ok=True)
        
        print("Testing BepiPred CLI help...")
        # Test 1: Check if BepiPred help works
        cmd_help = ["python", "bepipred3_CLI.py", "-h"]
        result_help = subprocess.run(cmd_help, capture_output=True, text=True, cwd=str(bepipred_path))
        
        if result_help.returncode == 0:
            print("✅ BepiPred help command works")
        else:
            print(f"❌ BepiPred help failed: {result_help.stderr}")
            return False
        
        print("Testing BepiPred prediction...")
        # Test 2: Run actual prediction
        cmd_pred = [
            "python", "bepipred3_CLI.py",
            "-i", fasta_file,
            "-o", os.path.abspath(output_dir),
            "-pred", "vt_pred"
        ]
        
        print(f"Command: {' '.join(cmd_pred)}")
        print(f"Working directory: {bepipred_path}")
        
        result_pred = subprocess.run(
            cmd_pred, 
            capture_output=True, 
            text=True, 
            cwd=str(bepipred_path),
            timeout=900  # 15 minutes for model download
        )
        
        print(f"Return code: {result_pred.returncode}")
        if result_pred.stdout:
            print(f"STDOUT: {result_pred.stdout}")
        if result_pred.stderr:
            print(f"STDERR: {result_pred.stderr}")
        
        if result_pred.returncode == 0:
            print("✅ BepiPred prediction completed")
            
            # Check output files
            output_files = os.listdir(output_dir)
            if output_files:
                print(f"✅ Output files created: {output_files}")
                
                # Look for results
                for file in output_files:
                    file_path = os.path.join(output_dir, file)
                    if os.path.getsize(file_path) > 0:
                        print(f"  - {file} ({os.path.getsize(file_path)} bytes)")
                        
                        # Show first few lines of text files
                        if file.endswith('.txt') or file.endswith('.csv'):
                            try:
                                with open(file_path, 'r') as f:
                                    lines = f.readlines()[:5]
                                    print(f"    Content preview: {lines}")
                            except:
                                pass
                
                return True
            else:
                print("⚠️  BepiPred ran but no output files found")
                return False
        else:
            print(f"❌ BepiPred prediction failed")
            return False
            
    except subprocess.TimeoutExpired:
        print("❌ BepiPred timed out")
        return False
    except Exception as e:
        print(f"❌ Error: {e}")
        return False
    finally:
        # Cleanup
        try:
            os.unlink(fasta_file)
        except:
            pass

if __name__ == "__main__":
    success = test_bepipred_direct()
    if success:
        print("\n🎉 BepiPred test successful!")
    else:
        print("\n💥 BepiPred test failed")
        print("\nTroubleshooting:")
        print("1. Make sure BepiPred dependencies are installed")
        print("2. Check if virtual environment needs to be activated")
        print("3. Try: cd tools/BepiPred3.0 && source venv/bin/activate && python bepipred3_CLI.py -h")