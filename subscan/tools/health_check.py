#!/usr/bin/env python3
"""
Health Check Tool for MutationScan Pipeline
Verifies all dependencies and system requirements.

Usage:
    python health_check.py              # Human-readable output
    python health_check.py --json       # JSON output for GUI integration
    python health_check.py --verbose    # Detailed diagnostic information
"""

import sys
import os
import subprocess
import platform
import shutil
import json
import argparse
from pathlib import Path
from typing import Dict, Any, Tuple


class Colors:
    """ANSI color codes for terminal output"""
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    RESET = '\033[0m'
    
    @staticmethod
    def strip():
        """Disable colors (for non-TTY output)"""
        Colors.GREEN = ''
        Colors.YELLOW = ''
        Colors.RED = ''
        Colors.BLUE = ''
        Colors.CYAN = ''
        Colors.BOLD = ''
        Colors.RESET = ''


def check_python_version() -> Dict[str, Any]:
    """Check Python version (requires 3.8+)"""
    version = sys.version_info
    version_str = f"{version.major}.{version.minor}.{version.micro}"
    meets_requirement = version >= (3, 8)
    
    return {
        'installed': True,
        'version': version_str,
        'meets_requirement': meets_requirement,
        'required': '3.8+',
        'status': 'good' if meets_requirement else 'error',
        'message': f"Python {version_str}" if meets_requirement else f"Python {version_str} (requires 3.8+)"
    }


def check_command(command: str, version_flag: str = '--version', timeout: int = 5) -> Dict[str, Any]:
    """Check if a command-line tool is available and get its version"""
    try:
        result = subprocess.run(
            [command, version_flag],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            text=True
        )
        
        # Extract version from output (first line usually)
        output = result.stdout + result.stderr
        version_line = output.split('\n')[0] if output else 'Unknown version'
        
        return {
            'installed': True,
            'version': version_line.strip(),
            'status': 'good',
            'message': f"{command} is available"
        }
    except FileNotFoundError:
        return {
            'installed': False,
            'version': None,
            'status': 'error',
            'message': f"{command} is not installed or not in PATH"
        }
    except subprocess.TimeoutExpired:
        return {
            'installed': False,
            'version': None,
            'status': 'error',
            'message': f"{command} command timed out"
        }
    except Exception as e:
        return {
            'installed': False,
            'version': None,
            'status': 'error',
            'message': f"Error checking {command}: {str(e)}"
        }


def check_emboss() -> Dict[str, Any]:
    """Check EMBOSS suite (specifically needle for alignments)"""
    # Check both needle and embossversion
    needle_result = check_command('needle', '-version')
    
    if needle_result['installed']:
        # Try to get EMBOSS version
        try:
            emboss_result = subprocess.run(
                ['embossversion'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=3,
                text=True
            )
            version = emboss_result.stdout.strip()
            needle_result['version'] = f"EMBOSS {version}"
        except:
            pass  # Keep needle version info
    
    return needle_result


def check_abricate() -> Dict[str, Any]:
    """Check ABRicate (AMR gene database tool)"""
    result = check_command('abricate', '--version')
    
    if result['installed']:
        # Check available databases
        try:
            db_result = subprocess.run(
                ['abricate', '--list'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=5,
                text=True
            )
            databases = [line.split()[0] for line in db_result.stdout.split('\n')[1:] if line.strip()]
            result['databases'] = databases
            result['message'] = f"ABRicate with {len(databases)} databases"
        except:
            pass
    
    return result


def check_disk_space(path: str = '.') -> Dict[str, Any]:
    """Check available disk space"""
    try:
        total, used, free = shutil.disk_usage(path)
        
        free_gb = free / (1024**3)
        total_gb = total / (1024**3)
        used_percent = (used / total) * 100
        
        # Warning if less than 10GB free
        status = 'good' if free_gb >= 10 else ('warning' if free_gb >= 5 else 'error')
        
        return {
            'free_gb': round(free_gb, 2),
            'total_gb': round(total_gb, 2),
            'used_percent': round(used_percent, 1),
            'status': status,
            'message': f"{round(free_gb, 1)}GB free" if status == 'good' else f"Low disk space: {round(free_gb, 1)}GB free"
        }
    except Exception as e:
        return {
            'free_gb': None,
            'total_gb': None,
            'used_percent': None,
            'status': 'error',
            'message': f"Error checking disk space: {str(e)}"
        }


def check_memory() -> Dict[str, Any]:
    """Check available memory (requires psutil if available)"""
    try:
        import psutil
        
        mem = psutil.virtual_memory()
        free_gb = mem.available / (1024**3)
        total_gb = mem.total / (1024**3)
        used_percent = mem.percent
        
        # Warning if less than 2GB free
        status = 'good' if free_gb >= 2 else ('warning' if free_gb >= 1 else 'error')
        
        return {
            'free_gb': round(free_gb, 2),
            'total_gb': round(total_gb, 2),
            'used_percent': round(used_percent, 1),
            'status': status,
            'message': f"{round(free_gb, 1)}GB free RAM" if status == 'good' else f"Low memory: {round(free_gb, 1)}GB free"
        }
    except ImportError:
        return {
            'free_gb': None,
            'total_gb': None,
            'used_percent': None,
            'status': 'warning',
            'message': "psutil not installed (optional)"
        }
    except Exception as e:
        return {
            'free_gb': None,
            'total_gb': None,
            'used_percent': None,
            'status': 'error',
            'message': f"Error checking memory: {str(e)}"
        }


def check_internet(test_url: str = 'https://www.ncbi.nlm.nih.gov', timeout: int = 5) -> Dict[str, Any]:
    """Check internet connectivity to NCBI"""
    try:
        import urllib.request
        
        urllib.request.urlopen(test_url, timeout=timeout)
        
        return {
            'connected': True,
            'test_url': test_url,
            'status': 'good',
            'message': 'Internet connection OK'
        }
    except Exception as e:
        return {
            'connected': False,
            'test_url': test_url,
            'status': 'warning',
            'message': f"Cannot reach {test_url}: {str(e)}"
        }


def check_python_packages() -> Dict[str, Any]:
    """Check essential Python packages"""
    packages = {}
    required = ['biopython', 'pandas', 'flask']
    optional = ['psutil', 'requests']
    
    for pkg in required + optional:
        try:
            if pkg == 'biopython':
                import Bio
                packages[pkg] = {'installed': True, 'version': Bio.__version__, 'required': True}
            elif pkg == 'pandas':
                import pandas
                packages[pkg] = {'installed': True, 'version': pandas.__version__, 'required': True}
            elif pkg == 'flask':
                import flask
                packages[pkg] = {'installed': True, 'version': flask.__version__, 'required': True}
            elif pkg == 'psutil':
                import psutil
                packages[pkg] = {'installed': True, 'version': psutil.__version__, 'required': False}
            elif pkg == 'requests':
                import requests
                packages[pkg] = {'installed': True, 'version': requests.__version__, 'required': False}
        except ImportError:
            packages[pkg] = {'installed': False, 'version': None, 'required': pkg in required}
    
    # Overall status
    missing_required = [pkg for pkg, info in packages.items() if info['required'] and not info['installed']]
    status = 'error' if missing_required else 'good'
    
    return {
        'packages': packages,
        'missing_required': missing_required,
        'status': status,
        'message': 'All required packages installed' if status == 'good' else f"Missing: {', '.join(missing_required)}"
    }


def run_health_checks() -> Dict[str, Any]:
    """Run all health checks and return results"""
    return {
        'system': {
            'platform': platform.system(),
            'platform_version': platform.version(),
            'machine': platform.machine(),
            'python_implementation': platform.python_implementation()
        },
        'python': check_python_version(),
        'packages': check_python_packages(),
        'emboss': check_emboss(),
        'abricate': check_abricate(),
        'disk': check_disk_space(),
        'memory': check_memory(),
        'internet': check_internet(),
        'overall_status': None  # Will be calculated
    }


def calculate_overall_status(results: Dict[str, Any]) -> str:
    """Calculate overall system health status"""
    # Check for any errors
    error_checks = [
        results['python']['status'] == 'error',
        results['packages']['status'] == 'error',
        results['emboss']['status'] == 'error',
        results['abricate']['status'] == 'error',
        results['disk']['status'] == 'error'
    ]
    
    if any(error_checks):
        return 'error'
    
    # Check for warnings
    warning_checks = [
        results['disk']['status'] == 'warning',
        results['memory']['status'] == 'warning',
        results['internet']['status'] == 'warning'
    ]
    
    if any(warning_checks):
        return 'warning'
    
    return 'good'


def print_human_readable(results: Dict[str, Any], verbose: bool = False):
    """Print health check results in human-readable format"""
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*70}{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.CYAN}MutationScan Pipeline - System Health Check{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*70}{Colors.RESET}\n")
    
    # System info
    if verbose:
        print(f"{Colors.BOLD}System Information:{Colors.RESET}")
        print(f"  Platform: {results['system']['platform']} ({results['system']['machine']})")
        print(f"  Python: {results['system']['python_implementation']}\n")
    
    # Python version
    py = results['python']
    status_icon = '✓' if py['status'] == 'good' else '✗'
    status_color = Colors.GREEN if py['status'] == 'good' else Colors.RED
    print(f"{status_color}{status_icon} Python:{Colors.RESET} {py['message']}")
    
    # Python packages
    pkg = results['packages']
    status_icon = '✓' if pkg['status'] == 'good' else '✗'
    status_color = Colors.GREEN if pkg['status'] == 'good' else Colors.RED
    print(f"{status_color}{status_icon} Python Packages:{Colors.RESET} {pkg['message']}")
    
    if verbose and pkg['packages']:
        for name, info in pkg['packages'].items():
            if info['installed']:
                req = " (required)" if info['required'] else " (optional)"
                print(f"    • {name} {info['version']}{req}")
            else:
                print(f"    • {name} {Colors.RED}NOT INSTALLED{Colors.RESET}")
    
    # EMBOSS
    emboss = results['emboss']
    status_icon = '✓' if emboss['status'] == 'good' else '✗'
    status_color = Colors.GREEN if emboss['status'] == 'good' else Colors.RED
    print(f"{status_color}{status_icon} EMBOSS:{Colors.RESET} {emboss['message']}")
    if verbose and emboss.get('version'):
        print(f"    Version: {emboss['version']}")
    
    # ABRicate
    abricate = results['abricate']
    status_icon = '✓' if abricate['status'] == 'good' else '✗'
    status_color = Colors.GREEN if abricate['status'] == 'good' else Colors.RED
    print(f"{status_color}{status_icon} ABRicate:{Colors.RESET} {abricate['message']}")
    if verbose and abricate.get('databases'):
        print(f"    Databases: {', '.join(abricate['databases'])}")
    
    # Disk space
    disk = results['disk']
    if disk['status'] == 'good':
        status_icon, status_color = '✓', Colors.GREEN
    elif disk['status'] == 'warning':
        status_icon, status_color = '⚠', Colors.YELLOW
    else:
        status_icon, status_color = '✗', Colors.RED
    print(f"{status_color}{status_icon} Disk Space:{Colors.RESET} {disk['message']}")
    if verbose and disk.get('total_gb'):
        print(f"    Total: {disk['total_gb']}GB, Used: {disk['used_percent']}%")
    
    # Memory
    mem = results['memory']
    if mem['status'] == 'good':
        status_icon, status_color = '✓', Colors.GREEN
    elif mem['status'] == 'warning':
        status_icon, status_color = '⚠', Colors.YELLOW
    else:
        status_icon, status_color = '✗', Colors.RED
    print(f"{status_color}{status_icon} Memory:{Colors.RESET} {mem['message']}")
    if verbose and mem.get('total_gb'):
        print(f"    Total: {mem['total_gb']}GB, Used: {mem['used_percent']}%")
    
    # Internet
    internet = results['internet']
    if internet['status'] == 'good':
        status_icon, status_color = '✓', Colors.GREEN
    elif internet['status'] == 'warning':
        status_icon, status_color = '⚠', Colors.YELLOW
    else:
        status_icon, status_color = '✗', Colors.RED
    print(f"{status_color}{status_icon} Internet:{Colors.RESET} {internet['message']}")
    
    # Overall status
    overall = results['overall_status']
    print(f"\n{Colors.BOLD}Overall Status:{Colors.RESET} ", end='')
    if overall == 'good':
        print(f"{Colors.GREEN}✓ System Ready{Colors.RESET}")
    elif overall == 'warning':
        print(f"{Colors.YELLOW}⚠ System Ready with Warnings{Colors.RESET}")
    else:
        print(f"{Colors.RED}✗ System Not Ready{Colors.RESET}")
    
    print(f"\n{Colors.CYAN}{'='*70}{Colors.RESET}\n")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='MutationScan Pipeline - System Health Check',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--json',
        action='store_true',
        help='Output results in JSON format (for GUI integration)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Show detailed diagnostic information'
    )
    
    args = parser.parse_args()
    
    # Disable colors for JSON output or non-TTY
    if args.json or not sys.stdout.isatty():
        Colors.strip()
    
    # Run health checks
    results = run_health_checks()
    results['overall_status'] = calculate_overall_status(results)
    
    # Output results
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print_human_readable(results, verbose=args.verbose)
    
    # Exit code based on overall status
    if results['overall_status'] == 'error':
        sys.exit(1)
    elif results['overall_status'] == 'warning':
        sys.exit(0)  # Warnings are OK to proceed
    else:
        sys.exit(0)


if __name__ == '__main__':
    main()
