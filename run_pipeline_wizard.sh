#!/bin/bash

# ==============================================================================
# GTEx Expression Analysis Pipeline Wizard
# ==============================================================================
# This script guides you through the gene expression analysis pipeline.
# It handles file paths, creates directories, and ensures a consistent workflow.
# ==============================================================================

# ANSI color codes for pretty output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Default values
DEFAULT_GTEX_TPM="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
DEFAULT_GTEX_ATTR="GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
BASE_OUTPUT_DIR="analysis_runs"

# Helper function to print headers
print_header() {
    echo -e "\n${BLUE}${BOLD}================================================================${NC}"
    echo -e "${BLUE}${BOLD}  $1${NC}"
    echo -e "${BLUE}${BOLD}================================================================${NC}\n"
}

# Helper function for user input with default
get_input() {
    local prompt="$1"
    local default="$2"
    local var_name="$3"
    local input
    
    if [ -z "$default" ]; then
        echo -e -n "${YELLOW}$prompt: ${NC}"
    else
        echo -e -n "${YELLOW}$prompt [${NC}${BOLD}$default${NC}${YELLOW}]: ${NC}"
    fi
    
    # Use read -e to enable Readline (provides filename completion)
    read -e input
    
    if [ -z "$input" ]; then
        input="$default"
    fi
    
    printf -v "$var_name" "%s" "$input"
}

# Helper function to confirm action
confirm() {
    local prompt="$1"
    echo -e -n "${YELLOW}$prompt [y/N]: ${NC}"
    read response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        return 0
    else
        return 1
    fi
}

# Helper to check if a file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo -e "${RED}Error: File '$1' not found.${NC}"
        return 1
    fi
    return 0
}

# ==============================================================================
# STEP 0: Welcome & Run Configuration
# ==============================================================================
clear
print_header "GTEx Expression Analysis Pipeline Wizard"
echo "This wizard will guide you through extracting, clustering, and analyzing"
echo "gene expression patterns from GTEx data."
echo ""

# Get Run Name
get_input "Enter a name for this analysis run (e.g., 'candidate_genes_v1')" "run_01" RUN_NAME

# Setup Directories
RUN_DIR="${BASE_OUTPUT_DIR}/${RUN_NAME}"
DATA_DIR="${RUN_DIR}/data"
RESULTS_DIR="${RUN_DIR}/results"
CONFIG_FILE="${RESULTS_DIR}/run_config.cfg"

echo -e "\n${GREEN}→ Setting up directories...${NC}"
mkdir -p "$DATA_DIR"
mkdir -p "$RESULTS_DIR"
echo "  Working Directory: $RUN_DIR"

# Check for existing configuration to allow resumption
CONFIG_LOADED=false
if [ -f "$CONFIG_FILE" ]; then
    print_header "Existing Project Found"
    echo -e "Found an existing configuration for project: ${BOLD}${RUN_NAME}${NC}"
    echo "This allows you to skip data extraction and use previous source file settings."
    confirm "Load previous project settings?"
    if [ $? -eq 0 ]; then
        source "$CONFIG_FILE"
        CONFIG_LOADED=true
        echo -e "${GREEN}✓ Settings loaded.${NC}"
    fi
fi

# ==============================================================================
# STEP 1: Data Extraction
# ==============================================================================
print_header "Step 1: Data Extraction"

# Define a local link for this specific run to allow easy resumption
RUN_EXPR_FILE="${DATA_DIR}/run_expression_matrix.csv"
SKIP_EXTRACT_INPUTS=false

# 1. Check for existing run data or loaded config to allow fast resume
if [ "$CONFIG_LOADED" = true ] && [ -f "$EXPRESSION_FILE" ]; then
    echo -e "${GREEN}✓ Using expression matrix from loaded config:${NC}"
    echo "  $EXPRESSION_FILE"
    SKIP_EXTRACT_INPUTS=true
elif [ -f "$RUN_EXPR_FILE" ]; then
    echo -e "${GREEN}✓ Found previously linked data for this run:${NC} $RUN_EXPR_FILE"
    confirm "Skip extraction inputs and use this existing data?"
    if [ $? -eq 0 ]; then
        EXPRESSION_FILE="$RUN_EXPR_FILE"
        SKIP_EXTRACT_INPUTS=true
    fi
fi

if [ "$SKIP_EXTRACT_INPUTS" = false ]; then
    # 2. Get Gene List
    while true; do
        get_input "Path to your gene list file (one gene symbol per line)" "candidate-genes.txt" GENES_FILE
        if check_file "$GENES_FILE"; then break; fi
    done

    # 3. Get GTEx Files (Needed to determine output filename version)
    get_input "Path to GTEx TPM file" "$DEFAULT_GTEX_TPM" GTEX_TPM_FILE
    if ! check_file "$GTEX_TPM_FILE"; then
        echo -e "${YELLOW}Warning: GTEx TPM file not found at default location. Please ensure it exists before running.${NC}"
    fi

    get_input "Path to GTEx Sample Attributes file" "$DEFAULT_GTEX_ATTR" GTEX_ATTR_FILE

    # 4. Determine GTEx Version for filename (e.g., v8, v10)
    # Matches 'v' followed by numbers in the filename
    GTEX_VERSION=$(echo "$GTEX_TPM_FILE" | grep -oE "v[0-9]+" | head -n 1)
    if [ -z "$GTEX_VERSION" ]; then
        GTEX_VERSION="custom"
    fi

    # 5. Construct shared expression filename
    GENES_BASENAME=$(basename "$GENES_FILE" .txt)
    SHARED_EXPRESSION_FILE="${BASE_OUTPUT_DIR}/${GENES_BASENAME}_${GTEX_VERSION}_expression_matrix.csv"

    echo -e "\nChecking for existing expression data for these genes and GTEx version (${GTEX_VERSION})..."
    
    DO_EXTRACT=true
    if [ -f "$SHARED_EXPRESSION_FILE" ]; then
        echo -e "${GREEN}✓ Found existing shared expression matrix:${NC} $SHARED_EXPRESSION_FILE"
        confirm "Do you want to reuse this existing data?"
        if [ $? -eq 0 ]; then
            DO_EXTRACT=false
        fi
    else
        echo "No pre-calculated expression matrix found for configuration: ${GENES_BASENAME} + ${GTEX_VERSION}"
    fi

    if [ "$DO_EXTRACT" = true ]; then
        echo -e "\n${GREEN}→ Running extraction...${NC}"
        echo "python extract_all_chunked.py --genes-file \"$GENES_FILE\" --output-file \"$SHARED_EXPRESSION_FILE\" --gtex-tpm-file \"$GTEX_TPM_FILE\" --sample-attributes \"$GTEX_ATTR_FILE\""
        
        python extract_all_chunked.py \
            --genes-file "$GENES_FILE" \
            --output-file "$SHARED_EXPRESSION_FILE" \
            --gtex-tpm-file "$GTEX_TPM_FILE" \
            --sample-attributes "$GTEX_ATTR_FILE"
            
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✓ Extraction complete.${NC}"
            echo -e "  Saved shared matrix to: $SHARED_EXPRESSION_FILE"
        else
            echo -e "${RED}❌ Extraction failed. Exiting.${NC}"
            exit 1
        fi
    else
        echo "Using existing shared expression matrix."
    fi
    
    # 6. Set the EXPRESSION_FILE for downstream use and Link it to the run directory
    EXPRESSION_FILE="$SHARED_EXPRESSION_FILE"
    
    # Create/Update the symlink so we can find it next time we run this wizard for this run_name
    # Use absolute path for link target to avoid relative path hell
    ABS_SHARED_PATH=$(readlink -f "$SHARED_EXPRESSION_FILE")
    ln -sf "$ABS_SHARED_PATH" "$RUN_EXPR_FILE"
    echo -e "  Linked to run directory: $RUN_EXPR_FILE"
fi

# ==============================================================================
# STEP 2: Clustering
# ==============================================================================
print_header "Step 2: Clustering Analysis"

confirm "Do you want to run/re-run clustering?"
if [ $? -eq 0 ]; then
    
    # Parameters - use loaded config values as defaults if available
    get_input "Cluster counts (k) to test (space-separated)" "${K_VALUES:-5 6 7 8 9}" K_VALUES
    
    echo -e "\n${BOLD}Normalization Options:${NC}"
    
    DEFAULT_Z_CONFIRM="n"
    [[ "$Z_FLAG" == "--z-score" ]] && DEFAULT_Z_CONFIRM="y"
    echo -e -n "${YELLOW}Apply Z-score standardization? (Recommended for tissue specificity) [${NC}${BOLD}${DEFAULT_Z_CONFIRM}${NC}${YELLOW}]: ${NC}"
    read response
    [[ -z "$response" ]] && response=$DEFAULT_Z_CONFIRM
    if [[ "$response" =~ ^[Yy]$ ]]; then Z_FLAG="--z-score"; else Z_FLAG=""; fi
    
    DEFAULT_SORT_CONFIRM="n"
    [[ "$SORT_FLAG" == "--cluster-sorted-heatmap" ]] && DEFAULT_SORT_CONFIRM="y"
    echo -e -n "${YELLOW}Use cluster-sorted heatmaps? (Easier to visualize clusters) [${NC}${BOLD}${DEFAULT_SORT_CONFIRM}${NC}${YELLOW}]: ${NC}"
    read response
    [[ -z "$response" ]] && response=$DEFAULT_SORT_CONFIRM
    if [[ "$response" =~ ^[Yy]$ ]]; then SORT_FLAG="--cluster-sorted-heatmap"; else SORT_FLAG=""; fi
    
    get_input "Statistic to use (median/mean)" "${STAT_METHOD:-median}" STAT_METHOD
    
    get_input "Minimum cluster size for zoomed heatmaps" "${MIN_CLUSTER_SIZE:-5}" MIN_CLUSTER_SIZE
    
    echo -e "\n${GREEN}→ Running clustering...${NC}"
    # Construct the command string for display
    CMD="python gtex_cluster_v2.py --expression-file \"$EXPRESSION_FILE\" --data-dir \"$RESULTS_DIR\" --k $K_VALUES --stat $STAT_METHOD $Z_FLAG $SORT_FLAG --min-cluster-size $MIN_CLUSTER_SIZE"
    echo "$CMD"
    
    python gtex_cluster_v2.py \
        --expression-file "$EXPRESSION_FILE" \
        --data-dir "$RESULTS_DIR" \
        --k $K_VALUES \
        --stat "$STAT_METHOD" \
        $Z_FLAG \
        $SORT_FLAG \
        --min-cluster-size "$MIN_CLUSTER_SIZE"
        
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Clustering complete.${NC}"
        echo -e "Outputs saved to: ${BLUE}$RESULTS_DIR${NC}"
    else
        echo -e "${RED}❌ Clustering failed.${NC}"
        # We allow continuing even if clustering fails, though subsequent steps might fail too.
        confirm "Continue anyway?" || exit 1
    fi
else
    echo "Skipping clustering step."
fi

# ==============================================================================
# STEP 2.5: Silhouette Visualization
# ==============================================================================
print_header "Step 2.5: Silhouette Visualization"

confirm "Generate Silhouette Plots (assess cluster quality)?"
if [ $? -eq 0 ]; then
    echo -e "\n${GREEN}→ Generating silhouette plots...${NC}"
    python visualise_silhouette_scores.py --data-dir "$RESULTS_DIR"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Plots generated in $RESULTS_DIR${NC}"
    else
        echo -e "${RED}❌ Visualization failed.${NC}"
    fi
fi

# ==============================================================================
# STEP 3: Enrichment Analysis
# ==============================================================================
print_header "Step 3: Enrichment Analysis"

# Check for cluster files
COUNT_CLUSTERS=$(find "$RESULTS_DIR" -name "*_genes.txt" | wc -l)
if [ "$COUNT_CLUSTERS" -eq 0 ]; then
    echo -e "${RED}No cluster gene files found in $RESULTS_DIR.${NC}"
    echo "Cannot run enrichment analysis without clusters."
    SKIP_ENRICHMENT=true
else
    echo -e "Found ${BOLD}$COUNT_CLUSTERS${NC} cluster gene lists."
    SKIP_ENRICHMENT=false
fi

if [ "$SKIP_ENRICHMENT" = false ]; then
    confirm "Run STRING Network Enrichment?"
    if [ $? -eq 0 ]; then
        echo -e "\n${GREEN}→ Running STRING enrichment...${NC}"
        python string_cluster_enrichment_report.py --data-dir "$RESULTS_DIR"
    fi

    confirm "Run Enrichr (Pathway/Ontology) Enrichment?"
    if [ $? -eq 0 ]; then
        echo -e "\n${GREEN}→ Running Enrichr analysis...${NC}"
        python enrichr_cluster_enrichment_report.py --data-dir "$RESULTS_DIR"
    fi
fi

# ==============================================================================
# SAVE CONFIGURATION
# ==============================================================================
# Save config before reporting so scripts can read it
CONFIG_FILE="${RESULTS_DIR}/run_config.cfg"
echo -e "${GREEN}→ Saving run configuration to ${CONFIG_FILE}...${NC}"

cat <<EOF > "$CONFIG_FILE"
# GTEx Pipeline Run Configuration
# Generated on $(date)

RUN_NAME="${RUN_NAME}"
GENES_FILE="$(readlink -f "$GENES_FILE")"
GTEX_TPM_FILE="$(readlink -f "$GTEX_TPM_FILE")"
GTEX_ATTR_FILE="$(readlink -f "$GTEX_ATTR_FILE")"
EXPRESSION_FILE="$(readlink -f "$EXPRESSION_FILE")"

# Clustering Parameters
K_VALUES="${K_VALUES}"
STAT_METHOD="${STAT_METHOD}"
Z_FLAG="${Z_FLAG}"
SORT_FLAG="${SORT_FLAG}"
MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE}"

# Paths
BASE_OUTPUT_DIR="$(readlink -f "$BASE_OUTPUT_DIR")"
RUN_DIR="$(readlink -f "$RUN_DIR")"
RESULTS_DIR="$(readlink -f "$RESULTS_DIR")"
DATA_DIR="$(readlink -f "$DATA_DIR")"
EOF

# ==============================================================================
# STEP 4: Reporting
# ==============================================================================
print_header "Step 4: Generating Reports"

if [ "$SKIP_ENRICHMENT" = false ]; then
    # Check if we actually have results to report on
    if [ -f "${RESULTS_DIR}/string_k_summary.csv" ] || [ -f "${RESULTS_DIR}/enrichr_k_summary.csv" ]; then
        echo -e "${GREEN}→ Generating Combined Enrichment Report...${NC}"
        python generate_combined_enrichment_report.py \
            --data-dir "$RESULTS_DIR" \
            --output-file "COMBINED_ENRICHMENT_REPORT.md"
    else
        echo -e "${YELLOW}Skipping Combined Report: No enrichment results found (did the API calls fail?)${NC}"
    fi
fi

if [ -f "$EXPRESSION_FILE" ] && [ "$COUNT_CLUSTERS" -gt 0 ]; then
    echo -e "${GREEN}→ Generating Biological Descriptions...${NC}"
    python generate_biological_descriptions.py \
        --data-dir "$RESULTS_DIR" \
        --expression-file "$EXPRESSION_FILE"
else
    echo "Skipping Biological Descriptions: Missing expression file or cluster files."
fi

# ==============================================================================
# FINISH
# ==============================================================================
print_header "Pipeline Complete!"

echo -e "All results are stored in: ${BLUE}${BOLD}$RUN_DIR${NC}"
echo ""
echo -e "${BOLD}Key Outputs:${NC}"
if [ -f "${RESULTS_DIR}/COMBINED_ENRICHMENT_REPORT.md" ]; then
    echo -e "  📄 Enrichment Report:   ${RESULTS_DIR}/COMBINED_ENRICHMENT_REPORT.md"
fi
if [ -f "${RESULTS_DIR}/CLUSTER_BIOLOGICAL_DESCRIPTIONS.md" ]; then
    echo -e "  📄 Biological Report:   ${RESULTS_DIR}/CLUSTER_BIOLOGICAL_DESCRIPTIONS.md"
fi
echo -e "  📊 Heatmaps & CSVs:     $RESULTS_DIR"

echo "\nYou can view the reports using any Markdown viewer or text editor."
