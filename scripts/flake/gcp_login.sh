#!/bin/bash

# GCP Login and Project Configuration Script
# This script handles GCP project setup and authentication

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if gcloud is installed
if ! command -v gcloud &> /dev/null; then
    print_error "gcloud CLI is not installed. Please install it first."
    exit 1
fi

# Check if we're already authenticated
if ! gcloud auth list --filter=status:ACTIVE --format="value(account)" | grep -q .; then
    print_status "No active GCP authentication found. Starting login process..."
    
    # Try to authenticate
    if gcloud auth login --no-launch-browser; then
        print_status "Successfully authenticated with GCP"
    else
        print_error "Failed to authenticate with GCP"
        exit 1
    fi
else
    print_status "Already authenticated with GCP"
fi

# Handle project configuration
if [ -f .gcp-project ]; then
    export GCP_PROJECT=$(cat .gcp-project)
    print_status "Using saved GCP_PROJECT=$GCP_PROJECT"
    
    # Check if we have access to this project
    if gcloud projects describe "$GCP_PROJECT" >/dev/null 2>&1; then
        gcloud config set project "$GCP_PROJECT"
        print_status "Successfully set GCP project to $GCP_PROJECT"
    else
        print_warning "Cannot access GCP project $GCP_PROJECT"
        print_warning "This may be due to insufficient authentication scopes on GCP VM"
        print_warning "The VM's service account may need additional scopes:"
        print_warning "  - https://www.googleapis.com/auth/cloud-platform"
        print_warning "  - https://www.googleapis.com/auth/secretmanager.readonly"
        print_warning "Attempting to continue with current configuration..."
    fi
else
    print_status "No saved GCP_PROJECT found."
    read -p "Enter your GCP project ID (or press Enter to skip): " entered_project
    
    if [ -n "$entered_project" ]; then
        export GCP_PROJECT=$entered_project
        echo "$GCP_PROJECT" > .gcp-project
        print_status "Saved GCP_PROJECT to .gcp-project"
        
        # Check if we have access to this project
        if gcloud projects describe "$GCP_PROJECT" >/dev/null 2>&1; then
            gcloud config set project "$GCP_PROJECT"
            print_status "Successfully set GCP project to $GCP_PROJECT"
        else
            print_warning "Cannot access GCP project $GCP_PROJECT"
            print_warning "This may be due to insufficient authentication scopes on GCP VM"
            print_warning "The VM's service account may need additional scopes:"
            print_warning "  - https://www.googleapis.com/auth/cloud-platform"
            print_warning "  - https://www.googleapis.com/auth/secretmanager.readonly"
            print_warning "Attempting to continue with current configuration..."
        fi
    else
        print_warning "No GCP project configured. Secret loading may fail."
        # Set a default empty project to prevent undefined variable issues
        export GCP_PROJECT=""
    fi
fi

print_status "GCP configuration complete!" 