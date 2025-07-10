#!/usr/bin/env bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <gcp-project> <secret-name>"
    exit 1
fi

GCP_PROJECT="$1"
SECRET_NAME="$2"

echo "[load-gcp-secret] Using GCP_PROJECT=$GCP_PROJECT"

# Check if GCP project is configured
if [ -z "$GCP_PROJECT" ]; then
    echo "[load-gcp-secret] WARNING: No GCP project configured. Skipping secret loading."
    exit 0
fi

# Check if gcloud is available
if ! command -v gcloud &> /dev/null; then
    echo "[load-gcp-secret] WARNING: gcloud CLI not found. Skipping secret loading."
    exit 0
fi

# Function to try fetching secret with current auth
try_fetch_secret() {
    local project="$1"
    local secret="$2"
    
    echo "[load-gcp-secret] Attempting to fetch secret '$secret' from project '$project'..."
    
    # Try to fetch the secret
    SECRET_VALUE=$(gcloud secrets versions access latest \
      --project "$project" \
      --secret "$secret" 2>&1)
    
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        echo "[load-gcp-secret] Successfully fetched secret!"
        return 0
    else
        echo "[load-gcp-secret] Failed to fetch secret: $SECRET_VALUE"
        return $exit_code
    fi
}

# Function to try fetching secret with Application Default Credentials
try_fetch_secret_with_adc() {
    local project="$1"
    local secret="$2"
    
    echo "[load-gcp-secret] Attempting to fetch secret with Application Default Credentials..."
    
    # Create a temporary gcloud configuration for ADC
    local temp_config="temp_adc_config_$$"
    
    # Create temporary config
    gcloud config configurations create "$temp_config" --no-activate >/dev/null 2>&1
    
    # Configure the temp config to use ADC
    gcloud config set auth/use_application_default_credentials true --configuration="$temp_config" >/dev/null 2>&1
    gcloud config set project "$project" --configuration="$temp_config" >/dev/null 2>&1
    
    # Try to fetch the secret using the temp config
    SECRET_VALUE=$(gcloud secrets versions access latest \
      --secret "$secret" \
      --configuration="$temp_config" 2>&1)
    
    local exit_code=$?
    
    # Clean up temporary config
    gcloud config configurations delete "$temp_config" --quiet >/dev/null 2>&1
    
    if [ $exit_code -eq 0 ]; then
        echo "[load-gcp-secret] Successfully fetched secret with ADC!"
        return 0
    else
        echo "[load-gcp-secret] Failed to fetch secret with ADC: $SECRET_VALUE"
        return $exit_code
    fi
}

# Check if we're authenticated with any method
AUTH_METHOD=""
if gcloud auth list --filter=status:ACTIVE --format="value(account)" | grep -q .; then
    AUTH_METHOD="service_account"
    echo "[load-gcp-secret] Using service account authentication"
elif gcloud auth application-default print-access-token >/dev/null 2>&1; then
    AUTH_METHOD="application_default"
    echo "[load-gcp-secret] Using application default credentials"
else
    echo "[load-gcp-secret] WARNING: No GCP authentication found. Skipping secret loading."
    exit 0
fi

# Try to fetch the secret
if try_fetch_secret "$GCP_PROJECT" "$SECRET_NAME"; then
    # Success! Parse and export the secret
    echo "$SECRET_VALUE" | jq -r 'to_entries[] | "export \(.key)=\(.value)"' > .env
    
    # Load into current shell
    set -a
    source .env
    set +a
    
    echo "[load-gcp-secret] Loaded $(jq length <<< "$SECRET_VALUE") keys into environment."
else
    # Check if it's an authentication scope issue
    if echo "$SECRET_VALUE" | grep -q "ACCESS_TOKEN_SCOPE_INSUFFICIENT"; then
        echo "[load-gcp-secret] This appears to be an authentication scope issue on GCP VM."
        echo "[load-gcp-secret] The VM's service account has insufficient scopes for Secret Manager."
        echo "[load-gcp-secret] Attempting to use application default credentials..."
        
        # Check if ADC is already set up
        if gcloud auth application-default print-access-token >/dev/null 2>&1; then
            echo "[load-gcp-secret] Application default credentials are already set up"
            if try_fetch_secret_with_adc "$GCP_PROJECT" "$SECRET_NAME"; then
                # Success with ADC! Parse and export the secret
                echo "$SECRET_VALUE" | jq -r 'to_entries[] | "export \(.key)=\(.value)"' > .env
                
                # Load into current shell
                set -a
                source .env
                set +a
                
                echo "[load-gcp-secret] Loaded $(jq length <<< "$SECRET_VALUE") keys into environment using ADC."
            else
                echo "[load-gcp-secret] Still failed to fetch secret with application default credentials"
                echo "[load-gcp-secret] Error: $SECRET_VALUE"
            fi
        else
            echo "[load-gcp-secret] Application default credentials not set up. Please run:"
            echo "[load-gcp-secret]   gcloud auth application-default login"
        fi
    else
        echo "[load-gcp-secret] Failed to fetch secret. Error: $SECRET_VALUE"
    fi
fi 