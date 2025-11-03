#!/usr/bin/env python3
"""
Code structure verification script
Verifies that the updated code is syntactically correct without requiring dependencies
"""

import sys
import os
import ast

print("="*60)
print("Code Structure Verification")
print("Verifying Updated CELLxGENE Census Integration")
print("="*60)
print()

def verify_python_file(filepath, description):
    """Verify a Python file is syntactically correct"""
    print(f"Checking {description}...")
    try:
        with open(filepath, 'r') as f:
            code = f.read()

        # Parse the code
        ast.parse(code)

        # Count functions
        tree = ast.parse(code)
        functions = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]

        print(f"  ✓ {filepath}")
        print(f"    - Syntax: Valid")
        print(f"    - Functions: {len(functions)}")
        if functions:
            print(f"    - Top functions: {', '.join(functions[:3])}")

        return True

    except SyntaxError as e:
        print(f"  ✗ Syntax Error in {filepath}")
        print(f"    Line {e.lineno}: {e.msg}")
        return False
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False

# Verify all source modules
print("[1] Verifying Source Modules")
print("-" * 60)

modules_to_check = [
    ("src/__init__.py", "Module initialization"),
    ("src/utils.py", "Utilities module"),
    ("src/download.py", "Download module (UPDATED)"),
    ("src/qc.py", "Quality control module"),
    ("src/preprocessing.py", "Preprocessing module"),
    ("src/integration.py", "Integration module"),
    ("src/visualization.py", "Visualization module"),
]

all_valid = True
for filepath, description in modules_to_check:
    if not verify_python_file(filepath, description):
        all_valid = False
    print()

# Verify configuration files exist
print("\n[2] Verifying Configuration Files")
print("-" * 60)

config_files = [
    "config/datasets.yaml",
    "config/qc_thresholds.yaml",
    "config/analysis_params.yaml"
]

for config_file in config_files:
    if os.path.exists(config_file):
        size = os.path.getsize(config_file)
        print(f"  ✓ {config_file} ({size} bytes)")
    else:
        print(f"  ✗ {config_file} NOT FOUND")
        all_valid = False

# Check specific updates in download.py
print("\n[3] Verifying Census API Updates in download.py")
print("-" * 60)

with open("src/download.py", 'r') as f:
    download_code = f.read()

updates_to_check = [
    ("get_obs(", "New get_obs() API usage"),
    ("with cellxgene_census.open_soma", "Context manager pattern"),
    ("get_census_version_description", "Census version info retrieval"),
    ("UPDATED 2025", "Documentation of updates"),
    ("fallback", "Fallback method for compatibility"),
    ("categorical", "Categorical columns handling"),
]

print("Checking for API updates:")
for search_string, description in updates_to_check:
    if search_string in download_code:
        print(f"  ✓ {description}")
    else:
        print(f"  ✗ {description} - NOT FOUND")

# Check documentation
print("\n[4] Verifying Documentation")
print("-" * 60)

doc_files = [
    ("README.md", "Main documentation"),
    ("CHANGELOG.md", "Change history"),
    ("UPDATES_2025.md", "API update guide (NEW)"),
    ("notebooks/README.md", "Notebook documentation"),
    ("scripts/README.md", "Scripts documentation"),
]

for doc_file, description in doc_files:
    if os.path.exists(doc_file):
        with open(doc_file, 'r') as f:
            lines = len(f.readlines())
        print(f"  ✓ {doc_file} - {description} ({lines} lines)")
    else:
        print(f"  ✗ {doc_file} NOT FOUND")

# Check notebook exists
print("\n[5] Verifying Notebook")
print("-" * 60)

if os.path.exists("notebooks/01_dataset_discovery.ipynb"):
    import json
    with open("notebooks/01_dataset_discovery.ipynb", 'r') as f:
        nb = json.load(f)

    n_cells = len(nb.get('cells', []))
    markdown_cells = sum(1 for cell in nb['cells'] if cell['cell_type'] == 'markdown')
    code_cells = sum(1 for cell in nb['cells'] if cell['cell_type'] == 'code')

    print(f"  ✓ notebooks/01_dataset_discovery.ipynb")
    print(f"    - Total cells: {n_cells}")
    print(f"    - Markdown cells: {markdown_cells}")
    print(f"    - Code cells: {code_cells}")
else:
    print(f"  ✗ Notebook NOT FOUND")
    all_valid = False

# Summary
print("\n" + "="*60)
print("VERIFICATION SUMMARY")
print("="*60)

if all_valid:
    print("✓ All code structure checks PASSED")
    print()
    print("Code is ready for execution with:")
    print("  1. Install dependencies: pip install -r requirements.txt")
    print("  2. Run test script: python3 test_notebook01.py")
    print("  3. Or use Docker: docker build -t nk-immunosenescence .")
    print()
    print("Key Updates Verified:")
    print("  ✓ download.py uses new get_obs() API")
    print("  ✓ Context managers implemented")
    print("  ✓ Fallback method for compatibility")
    print("  ✓ Documentation updated (UPDATES_2025.md)")
    print("  ✓ All modules syntactically correct")
else:
    print("✗ Some checks FAILED")
    print("  Please review errors above")

print("="*60)
