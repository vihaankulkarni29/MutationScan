#!/usr/bin/env python3
"""
Quick Final Test Runner
Runs essential validation tests for MutationScan pipeline
"""

import subprocess
import sys
from pathlib import Path

def run_quick_validation():
    """Run essential validation tests"""
    
    print("MutationScan Quick Validation Suite")
    print("="*50)
    
    tools_dir = Path(__file__).parent
    
    test_scripts = [
        ("Structure Check", "check_structure.py"),
        ("Dependencies", "check_dependencies.py"), 
        ("Smoke Tests", "smoke_test_dominos.py"),
        ("Full Domino Tests", "full_test_dominos.py"),
        ("Comprehensive Tests", "comprehensive_final_tests.py")
    ]
    
    results = {}
    
    for test_name, script in test_scripts:
        print(f"\n🔍 Running {test_name}...")
        
        try:
            result = subprocess.run(
                [sys.executable, script],
                cwd=tools_dir,
                capture_output=True,
                text=True
            )
            
            success = result.returncode == 0
            results[test_name] = success
            
            status = "✅ PASS" if success else "❌ FAIL"
            print(f"   {status}")
            
        except Exception as e:
            print(f"   ❌ ERROR: {e}")
            results[test_name] = False
    
    # Summary
    print(f"\n{'='*50}")
    print("VALIDATION SUMMARY")
    print(f"{'='*50}")
    
    passed = sum(results.values())
    total = len(results)
    
    for test_name, success in results.items():
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"  {test_name:20} | {status}")
    
    print(f"\nOverall: {passed}/{total} test suites passed")
    
    if passed == total:
        print("\n🚀 MutationScan pipeline is fully validated!")
        return True
    else:
        print("\n⚠️  Some validation tests failed. Check individual outputs.")
        return False

if __name__ == "__main__":
    success = run_quick_validation()
    sys.exit(0 if success else 1)