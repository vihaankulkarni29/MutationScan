#!/usr/bin/env python3
"""
MutationScan Web GUI Launcher
Starts the Flask web server and automatically opens the browser.

Usage:
    python subscan/gui/launch.py
    python subscan/gui/launch.py --port 8080
    python subscan/gui/launch.py --no-browser
    python subscan/gui/launch.py --debug
"""

import os
import sys
import time
import signal
import argparse
import webbrowser
import threading
from pathlib import Path

# Add parent directory to path to import app
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from subscan.gui.app import app


class Colors:
    """ANSI color codes for terminal output"""
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    END = '\033[0m'


def print_banner():
    """Print application banner"""
    banner = f"""
{Colors.BLUE}{Colors.BOLD}
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║              🧬  MutationScan Web Interface  🧬               ║
║                                                               ║
║         Comprehensive AMR Analysis Pipeline GUI               ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
{Colors.END}
"""
    print(banner)


def print_info(message):
    """Print info message"""
    print(f"{Colors.BLUE}ℹ{Colors.END}  {message}")


def print_success(message):
    """Print success message"""
    print(f"{Colors.GREEN}✓{Colors.END}  {message}")


def print_warning(message):
    """Print warning message"""
    print(f"{Colors.YELLOW}⚠{Colors.END}  {message}")


def print_error(message):
    """Print error message"""
    print(f"{Colors.RED}✗{Colors.END}  {message}")


def check_port_available(host, port):
    """Check if port is available"""
    import socket
    
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.settimeout(1)
    
    try:
        sock.bind((host, port))
        sock.close()
        return True
    except OSError:
        return False


def open_browser(url, delay=1.5):
    """Open browser after a short delay"""
    def _open():
        time.sleep(delay)
        print_info(f"Opening browser at {url}...")
        try:
            webbrowser.open(url)
            print_success("Browser opened successfully!")
        except Exception as e:
            print_warning(f"Could not open browser automatically: {e}")
            print_info(f"Please manually open: {url}")
    
    thread = threading.Thread(target=_open, daemon=True)
    thread.start()


def setup_signal_handlers(server_thread=None):
    """Setup graceful shutdown on Ctrl+C"""
    def signal_handler(sig, frame):
        print("\n")
        print_warning("Shutting down server...")
        print_info("Waiting for active connections to close...")
        time.sleep(1)
        print_success("Server stopped successfully!")
        print_info("Thank you for using MutationScan! 🧬")
        sys.exit(0)
    
    signal.signal(signal.SIGINT, signal_handler)
    if hasattr(signal, 'SIGTERM'):
        signal.signal(signal.SIGTERM, signal_handler)


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Launch MutationScan Web GUI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python subscan/gui/launch.py                    # Default: localhost:5000
  python subscan/gui/launch.py --port 8080        # Custom port
  python subscan/gui/launch.py --no-browser       # Don't open browser
  python subscan/gui/launch.py --debug            # Enable debug mode
  python subscan/gui/launch.py --host 0.0.0.0     # Listen on all interfaces
        """
    )
    
    parser.add_argument(
        '--port',
        type=int,
        default=5000,
        help='Port to run the server on (default: 5000)'
    )
    
    parser.add_argument(
        '--host',
        type=str,
        default='127.0.0.1',
        help='Host to bind to (default: 127.0.0.1). Use 0.0.0.0 for network access'
    )
    
    parser.add_argument(
        '--no-browser',
        action='store_true',
        help='Do not automatically open browser'
    )
    
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Run in debug mode (auto-reload on code changes)'
    )
    
    parser.add_argument(
        '--browser-delay',
        type=float,
        default=1.5,
        help='Delay before opening browser in seconds (default: 1.5)'
    )
    
    return parser.parse_args()


def validate_environment():
    """Validate that the environment is set up correctly"""
    print_info("Checking environment...")
    
    # Check Python version
    python_version = sys.version_info
    if python_version < (3, 8):
        print_error(f"Python 3.8+ required, found {python_version.major}.{python_version.minor}")
        return False
    print_success(f"Python {python_version.major}.{python_version.minor}.{python_version.micro}")
    
    # Check Flask is available
    try:
        import flask
        print_success(f"Flask {flask.__version__}")
    except ImportError:
        print_error("Flask not installed. Install with: pip install flask")
        return False
    
    # Check templates directory exists
    templates_dir = Path(__file__).parent / 'templates'
    if not templates_dir.exists():
        print_error(f"Templates directory not found: {templates_dir}")
        return False
    
    # Check required templates exist
    required_templates = ['setup.html', 'pipeline.html', 'results.html']
    missing_templates = []
    
    for template in required_templates:
        template_path = templates_dir / template
        if not template_path.exists():
            missing_templates.append(template)
    
    if missing_templates:
        print_error(f"Missing templates: {', '.join(missing_templates)}")
        return False
    
    print_success("All templates found")
    
    return True


def print_server_info(host, port, debug):
    """Print server information"""
    url = f"http://{host}:{port}"
    if host == '0.0.0.0':
        local_url = f"http://127.0.0.1:{port}"
        print(f"""
{Colors.GREEN}{Colors.BOLD}Server started successfully!{Colors.END}

{Colors.BOLD}Local Access:{Colors.END}
  {local_url}

{Colors.BOLD}Network Access:{Colors.END}
  Access from other devices on your network at:
  http://<your-ip-address>:{port}

{Colors.BOLD}Mode:{Colors.END} {'🐛 Debug (auto-reload enabled)' if debug else '🚀 Production'}

{Colors.BOLD}Controls:{Colors.END}
  • Press {Colors.BOLD}Ctrl+C{Colors.END} to stop the server
  • Server logs will appear below

{Colors.YELLOW}═══════════════════════════════════════════════════════════════{Colors.END}
""")
    else:
        print(f"""
{Colors.GREEN}{Colors.BOLD}Server started successfully!{Colors.END}

{Colors.BOLD}Access URL:{Colors.END}
  {url}

{Colors.BOLD}Mode:{Colors.END} {'🐛 Debug (auto-reload enabled)' if debug else '🚀 Production'}

{Colors.BOLD}Controls:{Colors.END}
  • Press {Colors.BOLD}Ctrl+C{Colors.END} to stop the server
  • Server logs will appear below

{Colors.YELLOW}═══════════════════════════════════════════════════════════════{Colors.END}
""")


def main():
    """Main entry point"""
    # Parse arguments
    args = parse_arguments()
    
    # Print banner
    print_banner()
    
    # Validate environment
    if not validate_environment():
        print_error("Environment validation failed. Please fix the issues above.")
        sys.exit(1)
    
    # Check port availability
    if not check_port_available(args.host, args.port):
        print_error(f"Port {args.port} is already in use!")
        print_info("Try a different port with --port flag, e.g.: --port 8080")
        sys.exit(1)
    
    print_success(f"Port {args.port} is available")
    
    # Setup signal handlers for graceful shutdown
    setup_signal_handlers()
    
    # Construct URL
    url = f"http://{'127.0.0.1' if args.host == '0.0.0.0' else args.host}:{args.port}"
    
    # Open browser (unless disabled)
    if not args.no_browser:
        open_browser(url, delay=args.browser_delay)
    else:
        print_info("Browser auto-open disabled. Manual access required.")
    
    # Print server information
    print_server_info(args.host, args.port, args.debug)
    
    # Start Flask server
    try:
        # Disable Flask's default startup message in production
        import logging
        if not args.debug:
            log = logging.getLogger('werkzeug')
            log.setLevel(logging.ERROR)
        
        # Run the app
        app.run(
            host=args.host,
            port=args.port,
            debug=args.debug,
            threaded=True,
            use_reloader=args.debug
        )
    except Exception as e:
        print_error(f"Failed to start server: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
