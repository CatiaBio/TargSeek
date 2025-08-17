#!/bin/bash

# download_blast_databases.sh
# Convenient script to download BLAST databases for epitope analysis
#
# Usage:
#   ./download_blast_databases.sh                    # Use config settings
#   ./download_blast_databases.sh --quick            # Essential databases only
#   ./download_blast_databases.sh --comprehensive    # All databases
#   ./download_blast_databases.sh --list             # Show available databases

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/config/config_download.yaml"

show_usage() {
    echo "Download BLAST databases for epitope cross-reactivity analysis"
    echo ""
    echo "Usage:"
    echo "  $0                      # Use current config settings"
    echo "  $0 --quick             # Essential databases only (human, viral, fungal, cow, mouse, ecoli)"
    echo "  $0 --comprehensive     # All available databases"
    echo "  $0 --list              # Show available databases and current settings"
    echo "  $0 --help              # Show this help"
    echo ""
    echo "Database categories:"
    echo "  • Standard: human, viral_proteins, fungal_proteins"
    echo "  • Livestock: cow, goat, horse, pig, chicken, buffalo, sheep"
    echo "  • Companion animals: dog, cat"
    echo "  • Model organisms: mouse, rat, rabbit"
    echo "  • Aquatic: salmon, zebrafish"
    echo "  • Plants: arabidopsis, rice, soybean"
    echo "  • Pathogens: ecoli, salmonella, staph"
    echo ""
    echo "Edit config/config_download.yaml to customize database selection"
}

show_databases() {
    echo "Available BLAST databases:"
    echo ""
    echo "Standard databases:"
    echo "  • human           - Human proteome (~25 MB)"
    echo "  • viral_proteins  - Viral proteins (~15 MB)"
    echo "  • fungal_proteins - Fungal proteins (~80 MB)"
    echo ""
    echo "Farm animals & livestock:"
    echo "  • cow             - Cow proteome (~22 MB)"
    echo "  • goat            - Goat proteome (~20 MB)"
    echo "  • horse           - Horse proteome (~20 MB)"
    echo "  • pig             - Pig proteome (~25 MB)"
    echo "  • chicken         - Chicken proteome (~17 MB)"
    echo "  • buffalo         - Water buffalo proteome (~20 MB)"
    echo "  • sheep           - Sheep proteome (~20 MB)"
    echo ""
    echo "Companion/domestic animals:"
    echo "  • dog             - Dog proteome (~25 MB)"
    echo "  • cat             - Cat proteome (~20 MB)"
    echo ""
    echo "Research model organisms:"
    echo "  • mouse           - Mouse proteome (~55 MB)"
    echo "  • rat             - Rat proteome (~30 MB)"
    echo "  • rabbit          - Rabbit proteome (~22 MB)"
    echo ""
    echo "Aquatic food organisms:"
    echo "  • salmon          - Atlantic salmon proteome (~35 MB)"
    echo "  • zebrafish       - Zebrafish proteome (~45 MB)"
    echo ""
    echo "Plant crops:"
    echo "  • arabidopsis     - Arabidopsis proteome (~35 MB)"
    echo "  • rice            - Rice proteome (~40 MB)"
    echo "  • soybean         - Soybean proteome (~56 MB)"
    echo ""
    echo "Bacterial pathogens:"
    echo "  • ecoli           - E. coli K-12 proteome (~1.3 MB)"
    echo "  • salmonella      - Salmonella proteome (~1.5 MB)"
    echo "  • staph           - Staphylococcus aureus proteome (~0.8 MB)"
    echo ""
    
    if [[ -f "$CONFIG_FILE" ]]; then
        echo "Current configuration in $CONFIG_FILE:"
        echo ""
        # Extract BLAST database settings
        python3 -c "
import yaml
try:
    with open('$CONFIG_FILE', 'r') as f:
        config = yaml.safe_load(f)
    
    blast_dbs = config.get('blast_databases', {})
    categories = {
        'Standard': ['human', 'viral_proteins', 'fungal_proteins'],
        'Livestock': ['cow', 'goat', 'horse', 'pig', 'chicken', 'buffalo', 'sheep'],
        'Companion': ['dog', 'cat'],
        'Model': ['mouse', 'rat', 'rabbit'],
        'Aquatic': ['salmon', 'zebrafish'],
        'Plants': ['arabidopsis', 'rice', 'soybean'],
        'Pathogens': ['ecoli', 'salmonella', 'staph']
    }
    
    for category, dbs in categories.items():
        enabled = [db for db in dbs if blast_dbs.get(db, False)]
        disabled = [db for db in dbs if not blast_dbs.get(db, False)]
        print(f'{category}:')
        if enabled:
            print(f'  ✅ Enabled: {", ".join(enabled)}')
        if disabled:
            print(f'  ❌ Disabled: {", ".join(disabled)}')
        print()

except Exception as e:
    print(f'Error reading config: {e}')
"
    fi
}

set_quick_config() {
    echo "Setting up quick configuration (essential databases only)..."
    
    # Backup original config
    cp "$CONFIG_FILE" "$CONFIG_FILE.backup.$(date +%Y%m%d_%H%M%S)"
    
    # Update config to enable only essential databases
    python3 -c "
import yaml

with open('$CONFIG_FILE', 'r') as f:
    config = yaml.safe_load(f)

# Set essential databases only
essential_dbs = ['human', 'viral_proteins', 'fungal_proteins', 'cow', 'mouse', 'ecoli']
blast_dbs = config.get('blast_databases', {})

for db in blast_dbs:
    blast_dbs[db] = db in essential_dbs

config['blast_databases'] = blast_dbs

with open('$CONFIG_FILE', 'w') as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False)
"
    echo "✅ Quick configuration applied"
    echo "Essential databases enabled: human, viral_proteins, fungal_proteins, cow, mouse, ecoli"
}

set_comprehensive_config() {
    echo "Setting up comprehensive configuration (all databases)..."
    
    # Backup original config
    cp "$CONFIG_FILE" "$CONFIG_FILE.backup.$(date +%Y%m%d_%H%M%S)"
    
    # Update config to enable all databases
    python3 -c "
import yaml

with open('$CONFIG_FILE', 'r') as f:
    config = yaml.safe_load(f)

# Enable all databases
blast_dbs = config.get('blast_databases', {})
for db in blast_dbs:
    blast_dbs[db] = True

config['blast_databases'] = blast_dbs

with open('$CONFIG_FILE', 'w') as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False)
"
    echo "✅ Comprehensive configuration applied"
    echo "All databases enabled"
}

run_download() {
    echo "Starting BLAST database download..."
    echo ""
    
    # Check if snakemake is available
    if ! command -v snakemake &> /dev/null; then
        echo "❌ Error: snakemake not found. Please install snakemake first:"
        echo "   conda install -c bioconda snakemake"
        exit 1
    fi
    
    # Check if makeblastdb is available
    if ! command -v makeblastdb &> /dev/null; then
        echo "❌ Error: makeblastdb not found. Please install BLAST+ tools:"
        echo "   conda install -c bioconda blast"
        exit 1
    fi
    
    # Run the download
    echo "Running: snakemake -s Snakefile_download all_blast_databases --cores 4"
    echo ""
    
    snakemake -s Snakefile_download all_blast_databases --cores 4
    
    if [[ $? -eq 0 ]]; then
        echo ""
        echo "🎉 BLAST database download completed successfully!"
        echo ""
        echo "Databases are located in: data/blast_databases/"
        echo "Download log: data/blast_databases/download.log"
        echo ""
        echo "💡 Usage example:"
        echo "   # Use in your epitope analysis scripts"
        echo "   blastp -query epitope.fasta -db data/blast_databases/human_proteome.fasta"
    else
        echo "❌ BLAST database download failed. Check the output above for errors."
        exit 1
    fi
}

# Main script logic
case "${1:-}" in
    --help|-h)
        show_usage
        ;;
    --list|-l)
        show_databases
        ;;
    --quick|-q)
        set_quick_config
        run_download
        ;;
    --comprehensive|-c)
        set_comprehensive_config
        run_download
        ;;
    "")
        run_download
        ;;
    *)
        echo "❌ Unknown option: $1"
        echo ""
        show_usage
        exit 1
        ;;
esac