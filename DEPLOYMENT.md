# Website Deployment Instructions

## Issue Identified

The website is **not deployed** because GitHub Pages has not been enabled in the repository settings. The GitHub Actions workflow builds successfully but fails during deployment with the following error:

```
Error: Failed to create deployment (status: 404)
Ensure GitHub Pages has been enabled: https://github.com/maxim-papusha/maxim-papusha.github.io/settings/pages
```

## Solution: Enable GitHub Pages

Follow these steps to enable GitHub Pages and deploy your website:

### Step 1: Enable GitHub Pages

1. Go to your repository on GitHub: https://github.com/maxim-papusha/maxim-papusha.github.io
2. Click on **Settings** (in the top navigation bar)
3. Scroll down and click on **Pages** (in the left sidebar under "Code and automation")
4. Under **Build and deployment**, find the **Source** section
5. Change the source from "Deploy from a branch" to **"GitHub Actions"**
6. The setting should now show: **Source: GitHub Actions**

### Step 2: Trigger a New Deployment

After enabling GitHub Pages with GitHub Actions as the source, you need to trigger a new deployment. You can do this in one of the following ways:

**Option A: Push a commit to main branch**
```bash
git commit --allow-empty -m "Trigger deployment"
git push origin main
```

**Option B: Re-run the failed workflow**
1. Go to the **Actions** tab in your repository
2. Find the failed "Deploy MkDocs to GitHub Pages" workflow run
3. Click on it
4. Click **Re-run all jobs**

### Step 3: Verify Deployment

1. Go back to **Settings → Pages**
2. After the workflow completes, you should see: **"Your site is live at https://maxim-papusha.github.io/"**
3. Click on the URL to view your deployed website

## Technical Details

### Repository Configuration

- **Site Generator**: MkDocs with Material theme
- **Build Command**: `mkdocs build --strict`
- **Output Directory**: `./site`
- **Deployment Method**: GitHub Actions workflow (`.github/workflows/deploy.yml`)

### Workflow Overview

The deployment workflow has two jobs:

1. **build**: Installs dependencies, builds the MkDocs site, and uploads it as an artifact
2. **deploy**: Deploys the artifact to GitHub Pages (only runs on pushes to main branch)

The workflow is configured with the necessary permissions:
```yaml
permissions:
  contents: read
  pages: write
  id-token: write
```

### What Was Verified

✅ Dependencies install correctly  
✅ MkDocs build completes successfully with `--strict` mode  
✅ Workflow configuration is correct  
✅ All necessary permissions are set  

❌ GitHub Pages is not enabled (requires manual configuration)

## After Deployment

Once GitHub Pages is enabled and the site is deployed, your blog will be accessible at:

**https://maxim-papusha.github.io/**

The site will automatically rebuild and redeploy every time you push changes to the main branch.

## Local Development

To preview the site locally before deploying:

```bash
# Install dependencies
pip install -r requirements.txt

# Serve the site locally (with live reload)
mkdocs serve

# Or build the site
mkdocs build --strict
```

The local server will be available at `http://127.0.0.1:8000/`
