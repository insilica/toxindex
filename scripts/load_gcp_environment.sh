#!/usr/bin/env bash

# Simple GCP secret loader
# Usage: ./load_gcp_environment.sh [secret-name]

# Embedded GCP project ID
GCP_PROJECT="873471276793"

if [ "$#" -gt 1 ]; then
    echo "Usage: $0 [secret-name]"
    echo "  If secret-name is omitted, loads all secrets from the project"
    exit 1
fi

SECRET_NAME="${1:-}"

# Clear the .env file at the start
> .env

if [ -z "$SECRET_NAME" ]; then
    echo "🔐 Loading all secrets from project '$GCP_PROJECT'..."
    
    # Get list of all secrets
    SECRETS=$(gcloud secrets list --project="$GCP_PROJECT" --format="value(name)")
    
    if [ -z "$SECRETS" ]; then
        echo "❌ No secrets found in project"
        exit 1
    fi
    
    # Count total secrets for progress
    TOTAL_SECRETS=$(echo "$SECRETS" | wc -l)
    CURRENT=0
    
    # Load each secret sequentially - GCP does not support multiple secrets being loaded at once
    for secret in $SECRETS; do
        CURRENT=$((CURRENT + 1))
        echo "📥 Loading secret ($CURRENT/$TOTAL_SECRETS): $secret"
        
        SECRET_VALUE=$(gcloud secrets versions access latest \
          --project "$GCP_PROJECT" \
          --secret "$secret" 2>/dev/null)
        
        if [ $? -eq 0 ] && [ -n "$SECRET_VALUE" ]; then
            # Check if it's JSON or single value
            if echo "$SECRET_VALUE" | jq . >/dev/null 2>&1; then
                # JSON secret - export each key/value
                echo "$SECRET_VALUE" | jq -r 'to_entries[] | "export \(.key)=\(.value)"' >> .env
                echo "✅ Loaded $(jq length <<< "$SECRET_VALUE") variables from $secret"
            else
                # Single value secret
                echo "export $secret=\"$SECRET_VALUE\"" >> .env
                echo "✅ Loaded single value: $secret"
            fi
        else
            echo "❌ Failed to load secret: $secret"
        fi
    done
    
    # Source the .env file after all secrets are loaded
    if [ -s .env ]; then
        set -a && source .env && set +a
        echo "🎉 Finished loading all secrets"
    else
        echo "⚠️  No secrets were loaded successfully"
    fi
else
    echo "🔐 Loading secret '$SECRET_NAME' from project '$GCP_PROJECT'..."
    
    # Fetch the secret
    SECRET_VALUE=$(gcloud secrets versions access latest \
      --project "$GCP_PROJECT" \
      --secret "$SECRET_NAME")
    
    if [ $? -eq 0 ]; then
        # Check if it's JSON or single value
        if echo "$SECRET_VALUE" | jq . >/dev/null 2>&1; then
            # JSON secret - export each key/value
            echo "$SECRET_VALUE" | jq -r 'to_entries[] | "export \(.key)=\(.value)"' > .env
            set -a && source .env && set +a
            echo "✅ Loaded $(jq length <<< "$SECRET_VALUE") environment variables"
        else
            # Single value secret
            echo "export $SECRET_NAME=\"$SECRET_VALUE\"" > .env
            export "$SECRET_NAME"="$SECRET_VALUE"
            echo "✅ Loaded single value: $SECRET_NAME"
        fi
    else
        echo "❌ Failed to load secret"
        exit 1
    fi
fi 