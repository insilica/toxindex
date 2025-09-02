#!/usr/bin/env bash
GCP_PROJECT="873471276793"
SECRET_NAME="OPENAI_API_KEY"

gcloud secrets versions access latest \
  --project "$GCP_PROJECT" \
  --secret "$SECRET_NAME"