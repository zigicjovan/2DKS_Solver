#!/bin/bash

# Universal GitHub upload script
# Works for both new and existing repos.
# Uploads only: data/, media/, and main_2DKS.m
# Usage:
#   ./gitpush.sh <project_path> <repo_url> "<commit_message>"
# Example:
#   ./gitpush.sh /home/user/2DKS https://github.com/username/2DKS.git "Update results"

PROJECT_DIR=$1
REPO_URL=$2
COMMIT_MSG=${3:-"Update from $(date)"}

if [ -z "$PROJECT_DIR" ]; then
  echo "❌ Usage: $0 <project_path> <repo_url> \"<commit_message>\""
  exit 1
fi

cd "$PROJECT_DIR" || { echo "❌ Project directory not found: $PROJECT_DIR"; exit 1; }

# Check if it's a git repo; if not, initialize
if [ ! -d ".git" ]; then
  echo "📦 Initializing new git repository..."
  git init
  git branch -M main
  if [ -z "$REPO_URL" ]; then
    echo "❌ Please provide a GitHub repo URL for first-time setup."
    exit 1
  fi
  git remote add origin "$REPO_URL"
fi

# Pull latest to prevent conflicts (ignore if remote not yet available)
git pull --rebase origin main 2>/dev/null

# Stage only specific directories and file
git add processing_functions/ solver_functions/ validation_scripts/ main_2DKS.m README.md 2>/dev/null

# Exit early if nothing changed
if git diff --cached --quiet; then
  echo "ℹ️  No changes detected in processing_functions/ solver_functions/ validation_scripts/ main_2DKS.m. README.md"
  exit 0
fi

# Commit and push
git commit -m "$COMMIT_MSG"
git push -u origin main

echo "✅ Successfully pushed selected files to GitHub!"

