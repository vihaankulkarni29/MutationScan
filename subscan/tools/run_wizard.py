#!/usr/bin/env python3
"""
MutationScan Interactive CLI Wizard
Guides users through pipeline configuration with interactive prompts.

Usage:
    python subscan/tools/run_wizard.py
"""

import os
import sys
import subprocess
from pathlib import Path
from datetime import datetime


class Colors:
    """ANSI color codes for terminal output"""
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    CYAN = '\033[96m'
    MAGENTA = '\033[95m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def print_banner():
    """Print wizard banner"""
    banner = f"""
{Colors.CYAN}{Colors.BOLD}
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                       ║
║            🧬  MutationScan Pipeline Setup Wizard  🧬                 ║
║                                                                       ║
║         Interactive Configuration for AMR Analysis Pipeline           ║
║                                                                       ║
╚═══════════════════════════════════════════════════════════════════════╝
{Colors.END}

{Colors.BOLD}Welcome!{Colors.END} This wizard will guide you through setting up your analysis.
Answer a few simple questions and we'll generate the command to run your pipeline.

"""
    print(banner)


def print_step(step_num, total_steps, title):
    """Print step header"""
    print(f"\n{Colors.BOLD}{Colors.BLUE}━━━ Step {step_num}/{total_steps}: {title} ━━━{Colors.END}\n")


def print_info(message):
    """Print info message"""
    print(f"{Colors.CYAN}ℹ{Colors.END}  {message}")


def print_success(message):
    """Print success message"""
    print(f"{Colors.GREEN}✓{Colors.END}  {message}")


def print_warning(message):
    """Print warning message"""
    print(f"{Colors.YELLOW}⚠{Colors.END}  {message}")


def print_error(message):
    """Print error message"""
    print(f"{Colors.RED}✗{Colors.END}  {message}")


def ask_question(question, default=None, validator=None):
    """Ask a question and get user input with optional validation"""
    while True:
        if default:
            prompt = f"{Colors.BOLD}{question}{Colors.END} [{Colors.CYAN}{default}{Colors.END}]: "
        else:
            prompt = f"{Colors.BOLD}{question}{Colors.END}: "
        
        answer = input(prompt).strip()
        
        # Use default if no answer provided
        if not answer and default:
            answer = default
        
        # Validate if validator provided
        if validator:
            valid, message = validator(answer)
            if not valid:
                print_error(message)
                continue
        
        return answer


def ask_yes_no(question, default=True):
    """Ask a yes/no question"""
    default_str = "Y/n" if default else "y/N"
    while True:
        answer = input(f"{Colors.BOLD}{question}{Colors.END} [{default_str}]: ").strip().lower()
        
        if not answer:
            return default
        
        if answer in ['y', 'yes']:
            return True
        elif answer in ['n', 'no']:
            return False
        else:
            print_error("Please answer 'y' or 'n'")


def ask_choice(question, choices, default=None):
    """Ask user to choose from a list of options"""
    print(f"\n{Colors.BOLD}{question}{Colors.END}")
    for idx, choice in enumerate(choices, 1):
        prefix = f"{Colors.GREEN}→{Colors.END}" if default and choice == default else " "
        print(f"  {prefix} {idx}. {choice}")
    
    while True:
        if default:
            prompt = f"\n{Colors.BOLD}Enter choice (1-{len(choices)}){Colors.END} [{Colors.CYAN}{choices.index(default)+1}{Colors.END}]: "
        else:
            prompt = f"\n{Colors.BOLD}Enter choice (1-{len(choices)}){Colors.END}: "
        
        answer = input(prompt).strip()
        
        # Use default if no answer
        if not answer and default:
            return default
        
        # Validate choice
        try:
            choice_idx = int(answer) - 1
            if 0 <= choice_idx < len(choices):
                return choices[choice_idx]
            else:
                print_error(f"Please enter a number between 1 and {len(choices)}")
        except ValueError:
            print_error("Please enter a valid number")


def validate_file(filepath):
    """Validate that a file exists"""
    if not filepath:
        return False, "File path cannot be empty"
    
    if not os.path.isfile(filepath):
        return False, f"File not found: {filepath}"
    
    return True, ""


def validate_email(email):
    """Validate email format"""
    if not email:
        return False, "Email cannot be empty"
    
    if '@' not in email:
        return False, "Invalid email format (must contain '@')"
    
    return True, ""


def validate_directory(dirpath):
    """Validate directory (create if doesn't exist)"""
    if not dirpath:
        return False, "Directory path cannot be empty"
    
    # Try to create if doesn't exist
    try:
        os.makedirs(dirpath, exist_ok=True)
        return True, ""
    except Exception as e:
        return False, f"Cannot create directory: {e}"


def create_file_from_input(filename):
    """Create a file from user input"""
    print(f"\n{Colors.YELLOW}Enter items (one per line). Press Ctrl+D (Linux/Mac) or Ctrl+Z then Enter (Windows) when done:{Colors.END}")
    
    items = []
    try:
        while True:
            line = input().strip()
            if line:
                items.append(line)
    except EOFError:
        pass
    
    if items:
        with open(filename, 'w') as f:
            for item in items:
                f.write(f"{item}\n")
        print_success(f"Created {filename} with {len(items)} items")
        return filename
    else:
        print_error("No items entered")
        return None


def suggest_common_genes():
    """Suggest common AMR genes"""
    suggestions = {
        'Fluoroquinolone Resistance': ['gyrA', 'parC', 'gyrB', 'parE'],
        'Beta-Lactamase': ['blaKPC', 'blaNDM', 'blaOXA', 'blaCTX', 'blaTEM', 'blaSHV'],
        'Aminoglycoside Resistance': ['aac', 'aph', 'ant'],
        'Tetracycline Resistance': ['tetA', 'tetB', 'tetC', 'tetM'],
        'Common Mix': ['gyrA', 'parC', 'blaKPC', 'blaNDM', 'aac', 'tetA']
    }
    
    print(f"\n{Colors.BOLD}Common gene categories:{Colors.END}")
    for idx, (category, genes) in enumerate(suggestions.items(), 1):
        print(f"  {idx}. {Colors.BOLD}{category}{Colors.END}: {', '.join(genes)}")
    
    return suggestions


def main():
    """Main wizard flow"""
    print_banner()
    
    # Configuration dictionary
    config = {}
    
    TOTAL_STEPS = 6
    
    # ========================================================================
    # STEP 1: Accession IDs
    # ========================================================================
    print_step(1, TOTAL_STEPS, "Genome Accession IDs")
    print_info("Provide NCBI accession IDs for the genomes you want to analyze")
    print_info("Example: NZ_CP107621, GCF_000005845.2")
    
    has_file = ask_yes_no("\nDo you have an accession list file?", default=True)
    
    if has_file:
        accessions_file = ask_question(
            "Enter path to accession list file",
            validator=validate_file
        )
        config['accessions'] = accessions_file
        
        # Show preview
        with open(accessions_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
            print_success(f"Loaded {len(lines)} accessions")
            print(f"  Preview: {', '.join(lines[:3])}{'...' if len(lines) > 3 else ''}")
    else:
        print_info("Let's create an accession list file")
        temp_file = f"accessions_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        accessions_file = create_file_from_input(temp_file)
        if accessions_file:
            config['accessions'] = accessions_file
        else:
            print_error("Cannot proceed without accessions. Exiting.")
            sys.exit(1)
    
    # ========================================================================
    # STEP 2: Target Genes
    # ========================================================================
    print_step(2, TOTAL_STEPS, "Target Genes")
    print_info("Specify which AMR genes you want to analyze")
    
    use_template = ask_yes_no("\nWould you like to use a common gene template?", default=True)
    
    if use_template:
        suggestions = suggest_common_genes()
        categories = list(suggestions.keys())
        chosen_category = ask_choice("\nSelect a category", categories, default=categories[0])
        
        # Create gene list file
        genes = suggestions[chosen_category]
        gene_file = f"genes_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        with open(gene_file, 'w') as f:
            for gene in genes:
                f.write(f"{gene}\n")
        
        print_success(f"Created gene list: {', '.join(genes)}")
        config['gene_list'] = gene_file
    else:
        has_gene_file = ask_yes_no("\nDo you have a gene list file?", default=True)
        
        if has_gene_file:
            gene_list_file = ask_question(
                "Enter path to gene list file",
                validator=validate_file
            )
            config['gene_list'] = gene_list_file
        else:
            print_info("Let's create a gene list file")
            temp_file = f"genes_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
            gene_list_file = create_file_from_input(temp_file)
            if gene_list_file:
                config['gene_list'] = gene_list_file
            else:
                print_error("Cannot proceed without gene list. Exiting.")
                sys.exit(1)
    
    # ========================================================================
    # STEP 3: Email for NCBI
    # ========================================================================
    print_step(3, TOTAL_STEPS, "Contact Information")
    print_info("NCBI requires an email address for API requests")
    
    config['email'] = ask_question(
        "Enter your email address",
        validator=validate_email
    )
    
    # ========================================================================
    # STEP 4: Species/Reference
    # ========================================================================
    print_step(4, TOTAL_STEPS, "Reference Sequences")
    print_info("Choose how to obtain wild-type reference sequences")
    
    common_species = [
        'Escherichia coli',
        'Klebsiella pneumoniae',
        'Pseudomonas aeruginosa',
        'Staphylococcus aureus',
        'Salmonella enterica',
        'Other (manual entry)'
    ]
    
    use_sepi = ask_yes_no("\nUse automatic reference download (recommended)?", default=True)
    
    if use_sepi:
        chosen_species = ask_choice("Select organism", common_species, default=common_species[0])
        
        if chosen_species == 'Other (manual entry)':
            config['sepi_species'] = ask_question("Enter species name (e.g., 'Acinetobacter baumannii')")
        else:
            config['sepi_species'] = chosen_species
        
        print_success(f"Will download references for: {config['sepi_species']}")
    else:
        ref_dir = ask_question(
            "Enter path to your reference sequence directory",
            validator=lambda p: (True, "") if os.path.isdir(p) else (False, "Directory not found")
        )
        config['user_reference_dir'] = ref_dir
    
    # ========================================================================
    # STEP 5: Database Selection
    # ========================================================================
    print_step(5, TOTAL_STEPS, "Database Source")
    print_info("Select which database to download genomes from")
    
    databases = ['NCBI', 'BV-BRC', 'EnteroBase', 'PATRIC']
    chosen_db = ask_choice("Select database", databases, default='NCBI')
    config['database'] = chosen_db.lower()
    
    # ========================================================================
    # STEP 6: Output & Advanced Options
    # ========================================================================
    print_step(6, TOTAL_STEPS, "Output Configuration")
    
    default_output = f"results/run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    config['output_dir'] = ask_question(
        "Enter output directory",
        default=default_output,
        validator=validate_directory
    )
    
    configure_advanced = ask_yes_no("\nConfigure advanced options (threads, etc.)?", default=False)
    
    if configure_advanced:
        config['threads'] = ask_question("Number of CPU threads", default="4")
        config['verbose'] = ask_yes_no("Enable verbose output?", default=False)
        config['open_report'] = ask_yes_no("Auto-open report when complete?", default=True)
    else:
        config['threads'] = "4"
        config['verbose'] = False
        config['open_report'] = True
    
    # ========================================================================
    # SUMMARY & COMMAND GENERATION
    # ========================================================================
    print(f"\n{Colors.BOLD}{Colors.GREEN}{'='*70}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.GREEN}Configuration Summary{Colors.END}")
    print(f"{Colors.BOLD}{Colors.GREEN}{'='*70}{Colors.END}\n")
    
    print(f"{Colors.BOLD}Accessions:{Colors.END} {config['accessions']}")
    print(f"{Colors.BOLD}Gene List:{Colors.END} {config['gene_list']}")
    print(f"{Colors.BOLD}Email:{Colors.END} {config['email']}")
    
    if 'sepi_species' in config:
        print(f"{Colors.BOLD}Species:{Colors.END} {config['sepi_species']} (auto-download)")
    else:
        print(f"{Colors.BOLD}Reference Dir:{Colors.END} {config['user_reference_dir']}")
    
    print(f"{Colors.BOLD}Database:{Colors.END} {config['database'].upper()}")
    print(f"{Colors.BOLD}Output Dir:{Colors.END} {config['output_dir']}")
    print(f"{Colors.BOLD}Threads:{Colors.END} {config['threads']}")
    
    # Generate command
    tools_dir = Path(__file__).parent
    cmd_parts = [
        sys.executable,
        str(tools_dir / 'run_pipeline.py'),
        '--accessions', config['accessions'],
        '--gene-list', config['gene_list'],
        '--email', config['email'],
        '--output-dir', config['output_dir'],
        '--threads', config['threads'],
    ]
    
    if 'sepi_species' in config:
        cmd_parts += ['--sepi-species', f'"{config["sepi_species"]}"']
    else:
        cmd_parts += ['--user-reference-dir', config['user_reference_dir']]
    
    if config.get('verbose'):
        cmd_parts.append('--verbose')
    
    if config.get('open_report'):
        cmd_parts.append('--open-report')
    
    command = ' '.join(cmd_parts)
    
    print(f"\n{Colors.BOLD}{Colors.CYAN}Generated Command:{Colors.END}")
    print(f"{Colors.YELLOW}{command}{Colors.END}\n")
    
    # Save command to file
    cmd_file = f"run_command_{datetime.now().strftime('%Y%m%d_%H%M%S')}.sh"
    with open(cmd_file, 'w') as f:
        f.write(f"#!/bin/bash\n")
        f.write(f"# MutationScan Pipeline Command\n")
        f.write(f"# Generated: {datetime.now().isoformat()}\n\n")
        f.write(command + "\n")
    
    print_success(f"Command saved to: {cmd_file}")
    
    # Ask to run now
    run_now = ask_yes_no(f"\n{Colors.BOLD}Run pipeline now?{Colors.END}", default=True)
    
    if run_now:
        print(f"\n{Colors.BOLD}{Colors.GREEN}Starting pipeline...{Colors.END}\n")
        print(f"{Colors.YELLOW}{'='*70}{Colors.END}\n")
        
        try:
            subprocess.run(cmd_parts, check=True)
            print(f"\n{Colors.YELLOW}{'='*70}{Colors.END}")
            print_success("Pipeline completed successfully!")
        except subprocess.CalledProcessError as e:
            print(f"\n{Colors.YELLOW}{'='*70}{Colors.END}")
            print_error(f"Pipeline failed with exit code {e.returncode}")
            sys.exit(1)
        except KeyboardInterrupt:
            print(f"\n{Colors.YELLOW}{'='*70}{Colors.END}")
            print_warning("Pipeline interrupted by user")
            sys.exit(130)
    else:
        print(f"\n{Colors.GREEN}Configuration complete!{Colors.END}")
        print(f"Run the pipeline later with:\n  {Colors.CYAN}bash {cmd_file}{Colors.END}\n")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print(f"\n\n{Colors.YELLOW}Wizard cancelled by user.{Colors.END}")
        sys.exit(130)
    except Exception as e:
        print_error(f"Unexpected error: {e}")
        sys.exit(1)
